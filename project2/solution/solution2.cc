/* 
Compiler Project 02: Auto-difference code generator from JSON to C++

Author: Bizhao Shi, Shaobo Zhou, Minggui Teng

Email: {shi_bizhao, telsazhou, minggui_teng}@pku.edu.cn

Date: June/19/2020

*/
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <assert.h>
#include <regex>
#include <set>
#include <stack>
#include <iterator>
#include <ctime>
#include "IR.h"
#include "IRMutator.h"
#include "IRVisitor.h"
#include "IRPrinter.h"
#include "type.h"
#include "json/json.h"
// #define DEBUG
using namespace Boost::Internal;

static int var_cnt = 0;
static int idx_cnt = 0;
static int imm_cnt = 0;
static int tmp_cnt = 0;
static int red_cnt = 0;
static int cid_cnt = 0;
const int typeNum = 4;

std::set<std::string> left_symbol;
std::set<std::string> both_symbol;
std::set<std::string> grad_set;
int global_data_type;
int max_step = 1;

class myIndex;
class myImmediate;
class myVariables;

const static std::string types_names[5] = {"float ", "int ", "unsigned int ", "std::string ", "handle "};

// A<16,32>[i+2, j+4]
class myVariables{
public:
	std::string name;				  //A
	int type;						  //float/int
	int dim;						  //[][]
	int is_left;					  //lef
	std::vector<unsigned long> shape;			  //1  A[i+j+3, i+j] a[][]
	std::vector<std::string> indexes; // index : idx_t
	std::set<std::string> symbols;	  //i+j+k i,j,k
	std::set<std::string> delta_symbols;
	int rw;
	bool need_grad;		//terminal
	bool contain_grad;	//non-terminal
	Expr expr;
	Expr grad_expr;
	myVariables(const myVariables & old);
	myVariables() {}
	myVariables(std::string s, int data_type, int flag); // flag is read or write sign
	void Print();
	void update_symbols(std::map<std::string, myIndex> &id);
	void Copy(const myVariables & old);
	bool check_pure(const myVariables & other);
};

class myIndex{
public:
	std::string name; //I+j+2
	unsigned long lb;
	unsigned long ub;
	int for_loop; // 0 :i+j j+2 2+j; 1: i  2: 2
	Expr dom;
	Expr slf;
	std::set<std::string> symbols;
	myIndex();
	int CheckIndex(std::string s);
	void Print();
	void IfCheck();
	bool operator<(const myIndex x) const{
		return (bool)name.compare(x.name); //从大到小排序
	}
};

class myImmediate{
public:
	std::string name; //imm_1/2/3
	int type;
	std::string value;
	myImmediate(){};
	void Print();
	Expr expr;
	Expr grad_expr;
};

static std::string pickone(std::string input, int &i){
	std::regex pattern("[a-zA-Z_][a-zA-Z0-9_]*");
	std::smatch result;
	std::string::const_iterator iterStart = input.begin() + i; //i+j+2
	std::string::const_iterator iterEnd = input.end();
	std::string temp;
	while (std::regex_search(iterStart, iterEnd, result, pattern)){
		temp = result[0];
		i = result[0].second - input.begin();
		break;
	}
	return temp;
}

static std::vector<std::string> inver(const std::string &input){
	std::vector<std::string> output;
	std::stack<char> operators;
	output.clear();
	for (int i = 0; i < (int)input.size(); i++){
		char ch = input[i];
		if (ch == '_' || (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')){
			std::string tmp = pickone(input, i);
			output.push_back(tmp);
			i--;
		}
		else if (ch == '+' || ch == '-'){
			while (!operators.empty() && (operators.top() == '*' || operators.top() == '/' || operators.top() == '-' || operators.top() == '+' || operators.top() == '%')){
				char tc[2] = {operators.top(), 0};
				output.push_back(std::string(tc));
				operators.pop();
			}
			operators.push(ch);
		}
		else if (ch == '*' || ch == '/' || ch == '%'){
			while (!operators.empty() && (operators.top() == '*' || operators.top() == '/' || operators.top() == '%')){
				char tc[2] = {operators.top(), 0};
				output.push_back(std::string(tc));
				operators.pop();
			}
			operators.push(ch);
			if (ch == '/' && input[i + 1] == '/'){
				i++;
			}
		}
		else if (ch == ')'){
			while (operators.top() != '(' && !operators.empty()){
				char tc[2] = {operators.top(), 0};
				output.push_back(std::string(tc));
				operators.pop();
			}
			operators.pop();
		}
		else if (ch == '('){
			operators.push(ch);
		}
	}
	while (!operators.empty()){
		char tc[2] = {operators.top(), 0};
		output.push_back(std::string(tc));
		operators.pop();
	}
	output.push_back(std::string("=\0"));
	return output;
}

static std::string pattern_replace(std::string prev, std::string curr, std::string src){
	std::regex pattern(prev);
	std::string result = regex_replace(src, pattern, curr);
	return result;
}

myVariables::myVariables(const myVariables & old){
	name.assign(old.name.begin(), old.name.end());				  //A
	type = old.type;						  //float/int
	dim = old.dim;						  //[][]
	is_left = old.is_left;					  //lef
	shape.assign(old.shape.begin(), old.shape.end());			  //1  A[i+j+3, i+j] a[][]
	indexes.assign(old.indexes.begin(), old.indexes.end()); // index : idx_t
	symbols.insert(old.symbols.begin(), old.symbols.end());	  //i+j+k i,j,k
	delta_symbols.insert(old.delta_symbols.begin(), old.delta_symbols.end());
	rw = old.rw;
	need_grad = old.need_grad;		//terminal
	contain_grad = old.contain_grad;	//non-terminal
}

void myVariables::Copy(const myVariables & old){
	name.assign(old.name.begin(), old.name.end());				  //A
	type = old.type;						  //float/int
	dim = old.dim;						  //[][]
	is_left = old.is_left;					  //lef
	shape.clear();
	shape.assign(old.shape.begin(), old.shape.end());			  //1  A[i+j+3, i+j] a[][]
	indexes.clear();
	indexes.assign(old.indexes.begin(), old.indexes.end()); // index : idx_t
	symbols.clear();
	symbols.insert(old.symbols.begin(), old.symbols.end());	  //i+j+k i,j,k
	delta_symbols.clear();
	delta_symbols.insert(old.delta_symbols.begin(), old.delta_symbols.end());
	rw = old.rw;
	need_grad = old.need_grad;		//terminal
	contain_grad = old.contain_grad;	//non-terminal
}

bool myVariables::check_pure(const myVariables & other){

	if(name.compare(other.name) != 0)	return 0;
	if(dim != other.dim)	return 0;
	for(int i = 0; i < dim; ++i){
		if(shape[i] != other.shape[i])	return 0;
		if(indexes[i].compare(other.indexes[i]) != 0)	return 0;
	}
	return 1;
}

myVariables::myVariables(std::string s, int data_type, int flag){
	name = s;
	type = data_type;
	rw = flag;
}

void myVariables::Print(){
	std::cout << "name: " << name << std::endl;
	std::cout << "type: " << types_names[type] << std::endl;
	std::cout << "read or write: " << rw << std::endl;
	std::cout << "dim: " << dim << std::endl;
		std::cout << "shape: ";
		for(int i = 0; i < (int)shape.size(); ++i){
			std::cout << "["<<shape[i] << "]";
		}
		for(int i = 0; i < (int)shape.size(); ++i){
			std::cout << "["<<shape[i] << "]";
		}
		std::cout << std::endl;
		std::cout << "symbols: ";
		for(auto & it : symbols){
			std::cout << it << ' ';
		}
		std::cout << std::endl;
		std::cout << "delta symbols: ";
		for(auto & it : delta_symbols){
			std::cout << it << ' ';
		}
		std::cout << std::endl;
}

void myVariables::update_symbols(std::map<std::string, myIndex> &id){
	if (indexes.size() == 0)
		return;
	for (int i = 0; i < (int)indexes.size(); ++i){
		std::regex pattern("[a-zA-Z_][a-zA-Z0-9_]*");
		std::smatch result;
		std::string::const_iterator iterStart = id[indexes[i]].name.begin(); //i+j+2
		std::string::const_iterator iterEnd = id[indexes[i]].name.end();
		std::string temp;
		while (std::regex_search(iterStart, iterEnd, result, pattern)){
			temp = result[0];
			symbols.insert(temp);
			iterStart = result[0].second;
		}
	}
}

myIndex::myIndex(){
	lb = 0;
	ub = 0x3f3f3f3f;
	symbols.clear();
}

int myIndex::CheckIndex(std::string s){
	int length = s.length();

	assert(length > 0); //is valid?

	if (s[0] == '_' || (s[0] >= 'a' && s[0] <= 'z') || (s[0] >= 'A' && s[0] <= 'Z')){ // i+2 / i
		for (int i = 1; i < length; ++i){
			if (!(s[i] == '_' || (s[i] >= 'a' && s[i] <= 'z') || (s[i] >= 'A' && s[i] <= 'Z') || (s[i] >= '0' && s[i] <= '9'))){
				return 0; // is a expression (i+2) or (i+j+2)
			}
		}
		return 1; // (i)
	}
	else if (s[0] >= '0' && s[0] <= '9'){
		for (int i = 1; i < length; ++i){
			if (!(s[0] >= '0' && s[0] <= '9')){
				return 0; // is a expression (2+j)
			}
		}
		return 2; //(2)
	}
	return -1;
}

void myIndex::Print(){
	std::cout << "name: " << name << std::endl;
	std::cout << "lower bound: " << lb << std::endl;
	std::cout << "upper bound: " << ub << std::endl;
	std::cout << "type: " << for_loop << std::endl
			  << std::endl;
}

void myImmediate::Print(){
	std::cout << "name: " << name << std::endl;
	std::cout << "type: " << types_names[type] << std::endl;
	std::cout << "value: " << value << std::endl
			  << std::endl;
}

static int judge_type(std::string s){
	for (int i = 0; i < typeNum; ++i){
		if (types_names[i].find(s) != std::string::npos){
			return i;
		}
	}
	std::cerr << "invalid data type!" << std::endl;
	assert(0);
}

static void process_io(const Json::Value &json_object, std::map<std::string, myVariables> &var,
					std::set<std::string>& grad, int start){
	global_data_type = judge_type(json_object["data_type"].asString());
	grad.clear();

		grad.insert(json_object["grad_to"][start].asString());

	//out
	int out_length = json_object["outs"].size();
	for (int i = 0; i < out_length; ++i){
		std::string tmp_name = json_object["outs"][i].asString();
		myVariables add_var = myVariables(tmp_name, global_data_type, -1);
		add_var.is_left = 1;
		var[tmp_name] = add_var; // read is 0
		var[tmp_name].need_grad = true;
		grad.insert(json_object["outs"][i].asString());
	}
	// in
	int in_length = json_object["ins"].size();
	for (int i = 0; i < in_length; ++i){
		std::string tmp_name = json_object["ins"][i].asString();
		myVariables add_var = myVariables(tmp_name, global_data_type, 1);
		add_var.is_left = -1;
		var[tmp_name] = add_var; // read is 0
		if(grad.find(tmp_name) != grad.end()){
			var[tmp_name].need_grad = true;
		}else{
			var[tmp_name].need_grad = false;
		}
	}
}

static std::vector<std::string> split_kernels(const Json::Value &json_object){
	std::string s = pattern_replace(" ", "", json_object["kernel"].asString());
	std::regex pattern(";");
	return std::vector<std::string>(std::sregex_token_iterator(s.begin(), s.end(), pattern, -1), std::sregex_token_iterator());
}

static bool CheckIsVariable(const std::map<std::string, myVariables> &var, std::string s){ // variables: A, B, C,  1A, 1B, 1C, ...
	auto it = var.find(s);
	if (it != var.end()){
		return 1;
	}
	else{
		return 0;
	}
}

static std::string CheckIsIndex(const std::map<std::string, myIndex> &indexes, std::string s){ // variables: A, B, C,  1A, 1B, 1C, ...
	int status = 1;
	for (auto &it : indexes){
		status = it.second.name.compare(s);
		if (status == 0){
			return it.first; //key
		}
	}
	std::string nf("NOT_FOUND");
	return nf;
}

struct para{
	bool flag;
	std::vector<unsigned long> shape;
};
static std::pair<int, std::string> process_variable(const std::string &kernel, int pos,
													std::map<std::string, myVariables> &var, std::map<std::string, myIndex> &index,
													int isleft){
	int tmp_pos = pos + 1; //A<16,32>[i,j] = (B[k,k]) + (B[i,k]*C[k,j])
	int klen = kernel.length();
	while (tmp_pos < klen){
		if (!(kernel[tmp_pos] == '_' ||
			  (kernel[tmp_pos] >= 'a' && kernel[tmp_pos] <= 'z') ||
			  (kernel[tmp_pos] >= 'A' && kernel[tmp_pos] <= 'Z') ||
			  (kernel[tmp_pos] >= '0' && kernel[tmp_pos] <= '9'))){
			break;
		}
		else{
			tmp_pos++;
		}
	}
	std::string tmp_name = kernel.substr(pos, tmp_pos - pos); // read until cannot read
	bool isVariable = CheckIsVariable(var, tmp_name);		  //A is variable
	if (isVariable){
		myVariables tmp_var;
		tmp_var.rw = 0;
		tmp_var.is_left = isleft;
		tmp_var.type = global_data_type;
		tmp_var.name = tmp_name;
		int start = kernel.find('<', tmp_pos - 1);
		int end = kernel.find('>', tmp_pos - 1);
		std::string tmp_range = kernel.substr(start, end - start + 1);
		// std::cout <<"     hhh" << ' ' << tmp_var.name << "  dim  "<<tmp_range <<std::endl;
		int range_length = tmp_range.length();
		int pos = 0;
		int number_length = 0;
		while (pos < range_length){
			while (tmp_range[pos] >= '0' && tmp_range[pos] <= '9'){
				number_length++;
				pos++;
			}
			if (number_length > 0){
				tmp_var.shape.push_back(atoi(tmp_range.substr(pos - number_length, number_length).c_str()));
				number_length = 0;
			}
			pos++;
		}
		tmp_var.dim = tmp_var.shape.size(); // dim is ok, shape is ok
		// std::cout <<"     hhh" << ' ' << tmp_var.name << ' ' << tmp_var.dim <<std::endl;
		// update shape
		if (var[tmp_var.name].dim == 0){
			var[tmp_var.name].dim = tmp_var.dim;
			for (int i = 0; i < tmp_var.dim; ++i){
				var[tmp_var.name].shape.push_back(tmp_var.shape[i]);
			}
		}

		if (tmp_var.dim == 1 && tmp_var.shape[0] == 1){															//alpha
			var_cnt++;												// expr cnt num
			std::string var_key = "var_" + std::to_string(var_cnt); //A[i] expr1[]

			tmp_var.indexes.push_back("");
			var[var_key] = tmp_var;
			// no symbol
			std::pair<int, std::string> tmp_pair(end + 1, var_key);
			return tmp_pair;
		}
		// [i, j]
		int last_end = end + 1;
		start = kernel.find('[', last_end);
		end = kernel.find(']', last_end);
		start++;
		end++;
		std::string tmp_index = kernel.substr(start, end - start - 1); //"A[i+j, j+2, i+k]"

		int index_length = tmp_index.length();
		pos = 0; //start_pos
		int end_pos = 0;
		for (int i = 0; i < tmp_var.dim; ++i){
			if (i < tmp_var.dim - 1){
				end_pos = tmp_index.find(',', pos);
			}
			else{
				end_pos = index_length;
			}
			if(tmp_var.dim == 1)	end_pos = index_length;
			std::string id_name = tmp_index.substr(pos, end_pos - pos);
			std::string root = CheckIsIndex(index, id_name);
			if (root.compare("NOT_FOUND") == 0){
				myIndex id;
				// check the valid
				id.name = id_name;
				id.ub = std::min(id.ub, tmp_var.shape[i]);
				id.lb = std::max(id.lb, (unsigned long)0);
				id.for_loop = id.CheckIndex(id.name);

				idx_cnt++; // expr cnt num
				std::string idx_key = "idx_" + std::to_string(idx_cnt);
				index[idx_key] = id; // add to the index map
				tmp_var.indexes.push_back(idx_key);
			}
			else{ //already exists
				index[root].ub = std::min(index[root].ub, tmp_var.shape[i]);
				index[root].lb = std::max(index[root].lb, (unsigned long)0);
				tmp_var.indexes.push_back(root); //idx_1/2/3/4
			}
			pos = end_pos + 1;
		}
		tmp_var.update_symbols(index);

		//
		if (isleft){
			for (auto &it : tmp_var.symbols){
				left_symbol.insert(it);
				both_symbol.insert(it);
			}
		}
		else{
			for (auto &it : tmp_var.symbols){ //tmp_var.symbols: i,j,   left_symbol: idx_1, idx_2
				both_symbol.insert(it);
				if (left_symbol.find(it) != left_symbol.end()){
					continue;
				}
				else{
					tmp_var.delta_symbols.insert(it);
				}
			}
		}
		var_cnt++;												// expr cnt num
		std::string var_key = "var_" + std::to_string(var_cnt); //A[i] expr1[]
		tmp_var.need_grad = false;
		// if(grad_set.find(tmp_var.name) != grad_set.end())	tmp_var.need_grad = true;
		// else
		// {
		// 	tmp_var.need_grad = false;
		// }
		
		var[var_key] = tmp_var;
		std::pair<int, std::string> tmp_pair(end, var_key);
		return tmp_pair;
	}
	else{ //index
	}
	std::pair<int, std::string> tmp_pair(-1, "");
	return tmp_pair;
}

static std::pair<int, std::string> process_immediate(const std::string &kernel,
													 int cur_pos, std::map<std::string, myImmediate> &imm){
	myImmediate number;
	number.type = global_data_type;
	int rv = cur_pos;
	switch (global_data_type){
	case 0:{ //float
		std::smatch result;
		std::regex pattern("\\d+(\\.\\d*)?|\\.\\d+"); //匹配四个数字
		std::string::const_iterator iterStart = kernel.begin() + cur_pos;
		std::string::const_iterator iterEnd = kernel.end();
		std::string temp;
		while (std::regex_search(iterStart, iterEnd, result, pattern)){
			rv += result[0].second - iterStart;
			number.value = result[0];
			break;
		}
		break;
	}
	case 1:{ // int
		std::smatch result;
		std::regex pattern("\\d+(\\.\\d*)?|\\.\\d+"); //匹配四个数字
		std::string::const_iterator iterStart = kernel.begin() + cur_pos;
		std::string::const_iterator iterEnd = kernel.end();
		std::string temp;
		while (std::regex_search(iterStart, iterEnd, result, pattern)){
			rv += result[0].second - iterStart;
			std::string tmp = result[0];
			number.value = std::to_string(atoi(tmp.c_str()));
			break;
		}
		break;
	}
	case 2:{ // unsigned int
		unsigned int uvalue;
		sscanf(kernel.substr(cur_pos).c_str(), "%u", &uvalue);
		number.value = std::to_string(uvalue);
		rv += number.value.length();
		break;
	}
	case 3:{ //string
		char cvalue[1024];
		sscanf(kernel.substr(cur_pos).c_str(), "%s", cvalue);
		number.value = cvalue;
		rv += number.value.length();
		break;
	}
	default:
		break;
		assert(0);
	};

	imm_cnt++;
	std::string imm_key = "imm_" + std::to_string(imm_cnt);
	number.name = imm_key;
	imm[imm_key] = number;
	std::pair<int, std::string> tmp_pair(rv, imm_key);
	return tmp_pair;
}
// "A[i][j] = ,..."
// "var_1 = var_2 * var_3 + imm_1 * var_4"

static std::string process_kernel(const std::string &kernel,
								  std::map<std::string, myVariables> &var,
								  std::map<std::string, myIndex> &index,
								  std::map<std::string, myImmediate> &imm){
	std::string my_ir;
	int length = kernel.length();
	int cur_pos = 0;
	int isleft = 1;
	while (cur_pos < length){ // suffix expr
		char cur_char = kernel[cur_pos];
		if (cur_char == '_' || (cur_char >= 'a' && cur_char <= 'z') || (cur_char >= 'A' && cur_char <= 'Z')){																								//variable state
			std::pair<int, std::string> mypair = process_variable(kernel, cur_pos, var, index, isleft); //after processing, the cur position updates
			cur_pos = mypair.first;
			my_ir += mypair.second; // var_1
		}
		else if (cur_char == '.' || (cur_char >= '0' && cur_char <= '9')){
			std::pair<int, std::string> mypair = process_immediate(kernel, cur_pos, imm);
			cur_pos = mypair.first;
			my_ir += mypair.second;
		}
		else if (cur_char == '='){
			isleft = 0;
			my_ir += cur_char;
			cur_pos++;
		}
		else{
			my_ir += cur_char;
			cur_pos++;
		}
	}

	return my_ir;
}

static void create_expr(std::map<std::string, myVariables> &var, std::map<std::string, myImmediate> &imm,
						std::map<std::string, myIndex> &id){
	Type index_type = Type::int_scalar(32);
	for (auto &it : id){
		it.second.dom = Dom::make(index_type, (int)it.second.lb, (int)it.second.ub);
		it.second.slf = Index::make(index_type, it.second.name, it.second.dom, IndexType::Spatial);
	}

	Type data_type; //giao
	Expr tmp_zero_expr;
	switch (global_data_type){
	case 0: // float
		data_type = Type::float_scalar(32);
		tmp_zero_expr = FloatImm::make(data_type, 0);
		break;
	case 1:
		data_type = Type::int_scalar(32);
		tmp_zero_expr = IntImm::make(data_type, 0);
		break;
	case 2:
		data_type = Type::uint_scalar(32);
		tmp_zero_expr = UIntImm::make(data_type, 0);
		break;
	default:
		break;
	}
	for (auto &it : var){
		if (it.first[0] != 'v'){
			continue;
		}	
		if (it.second.dim == 1 && it.second.shape[0] == 1){
			it.second.expr = Var::make(data_type, it.second.name, {}, {1});
			if(it.second.need_grad == true || it.second.is_left == true){
				it.second.grad_expr = Var::make(data_type, "d"+it.second.name, {},{1});
			}else{
				it.second.grad_expr = tmp_zero_expr;
			}
		}
		else{
			std::vector<Expr> tmp_id;
			std::vector<long unsigned int> my_shape;
			for (int dim_i = 0; dim_i < it.second.dim; ++dim_i){
				tmp_id.push_back(id[it.second.indexes[dim_i]].slf);
				my_shape.push_back((long unsigned int)(it.second.shape[dim_i]));
			}
			it.second.expr = Var::make(data_type, it.second.name, tmp_id, my_shape);
			if(it.second.need_grad == true || it.second.is_left == true){
				it.second.grad_expr = Var::make(data_type, "d"+it.second.name, tmp_id, my_shape);
			}else{
				it.second.grad_expr = tmp_zero_expr;
			}
		}
	}
	for (auto &it : imm){
		switch (it.second.type){
		case 0:{ //float
			Type imm_type = Type::float_scalar(32);
			it.second.expr = FloatImm::make(imm_type, atof(it.second.value.c_str()));
			it.second.grad_expr = tmp_zero_expr;
			break;
		}
		case 1:{
			Type imm_type = Type::int_scalar(32);
			it.second.expr = IntImm::make(imm_type, atoi(it.second.value.c_str()));
			it.second.grad_expr = tmp_zero_expr;
			break;
		}
		case 2:{
			Type imm_type = Type::uint_scalar(32);
			it.second.expr = UIntImm::make(imm_type, (unsigned int)atoi(it.second.value.c_str()));
			it.second.grad_expr = tmp_zero_expr;
			break;
		}
		case 3:{
			// Type imm_type = Type::string_scalar(32);
			// Expr int_expr = UIntImm::make(imm_type, (unsigned int)atoi(pointer->value.c_str()));
			// s.push(int_expr);
		}
		default:{
			std::cerr << "invalid data type detected!" << std::endl;
			break;
		}
		}
	}
}

static Expr index_based_cond_gen(std::vector<myIndex> &vec, int n){

	Type index_type = Type::int_scalar(32); // not giao!

	Expr lee = Compare::make(index_type, CompareOpType::GE, vec[n].slf, IntImm::make(index_type, vec[n].lb));
	Expr gee = Compare::make(index_type, CompareOpType::LT, vec[n].slf, IntImm::make(index_type, vec[n].ub));
	Expr bee = Binary::make(index_type, BinaryOpType::And, lee, gee);

	int size = vec.size();
	if (n == size - 1)
		return Binary::make(index_type, BinaryOpType::And, lee, gee);
	else
		return Binary::make(index_type, BinaryOpType::And, bee, index_based_cond_gen(vec, n + 1));
}

static Expr expr_gen(Expr &in_expr1, Expr &in_expr2, int opcode){
	Type data_type; //giao
	switch (global_data_type){
	case 0: // float
		data_type = Type::float_scalar(32);
		break;
	case 1:
		data_type = Type::int_scalar(32);
		break;
	case 2:
		data_type = Type::uint_scalar(32);
		break;
	default:
		break;
	}
	switch (opcode){
	case 0: //mult
		return Binary::make(data_type, BinaryOpType::Mul, in_expr1, in_expr2);
	case 1: //div
		return Binary::make(data_type, BinaryOpType::Div, in_expr1, in_expr2);
	case 2: //mod
		return Binary::make(data_type, BinaryOpType::Mod, in_expr1, in_expr2);
	}
}

static Expr pass_sub_str(std::vector<std::string> &str, //suffix
						 std::map<std::string, myIndex> &id,
						 std::map<std::string, myVariables> &var,
						 std::map<std::string, myImmediate> &imm,
						 Type data_type, Expr tmp_zero_expr, bool* grad_flag){
	
	Expr all = tmp_zero_expr;
	std::string op_string = "*/%";
	std::stack<Expr> s;

	// update the max step
	int current_step = 0;
	for(auto & it : str){
		if(it[0] == 'v'){
			current_step += (var[it].need_grad == true);
		}
	}
	max_step = std::max(current_step, max_step);
	// update the grad_flag
	bool cur_grad_flag = false;
	for(auto & it : str){
		if(it[0] == 'v'){
			cur_grad_flag |= var[it].need_grad;
		}
		if(it[0] == 'r'){
			cur_grad_flag |= var[it].contain_grad;
		}
	}
	*grad_flag = cur_grad_flag;

	max_step = std::max(current_step, max_step);
	for(int i = 0; i < str.size(); ++i){
		// std::cout << "outer loop: " << i << std::endl;
		if (op_string.find(str[i]) != std::string::npos)	continue;
		if (str[i].find("imm_") != std::string::npos)	continue;
		if (str[i].find("red_") != std::string::npos && var[str[i]].contain_grad == false)	continue;
		if(str[i].find("var_") != std::string::npos && var[str[i]].need_grad == false)	continue;
		std::string target = str[i];	//grad object
		// std::cout << "this loop grad index: " << str[i] << std::endl;
		for (auto &it : str){ //var_1 *
			// std::cout <<"inter loop: " << it <<std::endl;
			if (op_string.find(it) != std::string::npos){
				Expr top1 = s.top();
				s.pop();
				Expr top2 = s.top();
				s.pop();
				// std::cout <<"expr gen before" << std::endl;
				s.push(expr_gen(top2, top1, op_string.find(it)));
				// std::cout <<"expr gen end" << std::endl;
			}
			else{
				if(it.compare(target) == 0){
					if (it.find("var_") != std::string::npos){
						// std::cout <<"case10hhh"<<std::endl;
						s.push(var[it].grad_expr);
					}
					else if (it.find("red_") != std::string::npos){
						s.push(var[it].expr);
					}
					else if (it.find("imm_") != std::string::npos){
						s.push(tmp_zero_expr);
					}	
				}	
				else{
					if (it.find("var_") != std::string::npos){
						// if(var[it].need_grad == true)
						// 	s.push(var[it].grad_expr);
						// else{
							s.push(var[it].expr);
	//						}
						}
					else if (it.find("red_") != std::string::npos){
						s.push(var[it].expr);
					}
					else if (it.find("imm_") != std::string::npos){
						s.push(imm[it].expr);
					}
				}
			}
		}
		Expr top = s.top();
		s.pop();
		all = Binary::make(data_type, BinaryOpType::Add, all, top);
	}
	return all;
}

static Expr sub_reduction(std::vector<std::string> &str, //suffix
						  std::map<std::string, myVariables> &var,
						  std::map<std::string, myIndex> &id,
						  std::map<std::string, myImmediate> &imm,
						  std::vector<Stmt> &core_code,
						  Type data_type, Expr tmp_zero_expr, bool* grad_flag){
	// red_1
	if (str.size() == 1 && str[0].find("red_") != std::string::npos){
		return var[str[0]].expr;
	}
	else{ // var_1 imm_
		// std::cout << "before expr " << std::endl;
		Expr core_expr = pass_sub_str(str, id, var, imm, data_type, tmp_zero_expr, grad_flag);
	// std::cout << "after expr " << std::endl;
		// tmp_var
		std::set<std::string> loop_se;
		loop_se.clear();
		int str2_size = str.size();
		for (int i = 0; i < str2_size; ++i){ // var1 imm1 * var2 *
			if (str[i].find("var_") != std::string::npos){
				for (auto &it : var[str[i]].delta_symbols){ //idx_1/2/3
					loop_se.insert(it);
				}
			}
		}

		// if (loop_se.size() == 0)
		// 	return core_expr;

		tmp_cnt++;
		std::string tmp_var_name = "tmp_";
		tmp_var_name += std::to_string(tmp_cnt);
		Expr expr_tmp = Var::make(data_type, tmp_var_name, {}, {1});									   // tmp_x
		Expr expr_tmp_first = Var::make(data_type, types_names[global_data_type] + tmp_var_name, {}, {1}); // float tmp_x;
		Stmt tmp_init = Move::make(expr_tmp_first, IntImm::make(data_type, 0), MoveType::MemToMem);
		core_code.push_back(tmp_init);
		std::vector<Expr> lpv;

		lpv.clear();

		for (auto &it : loop_se){
			for (auto &itt : id){
				if (itt.second.name.compare(it) == 0){
					lpv.push_back(itt.second.slf);
				}
			}
		}

		std::set<myIndex> idx_se;
		idx_se.clear();
		int str_size = str.size();
		for (int i = 0; i < str_size; ++i){ // var1 imm1 * var2 *
			if (str[i].find("var_") != std::string::npos){
				std::cout << str[i] << std::endl;
				if (var[str[i]].dim == 1 && var[str[i]].shape[0] == 1)
					continue;
				for (auto &it : var[str[i]].indexes){ //idx_1/2/3
				std::cout << id[it].name << std::endl;
					if (id[it].for_loop > 0)
						continue;
					// std::cout <<"insert id:                          " << id[it].name << std::endl;
					idx_se.insert(id[it]);
				}
			}
		}
		std::vector<myIndex> idv(idx_se.begin(), idx_se.end());

		Expr eadd = Binary::make(data_type, BinaryOpType::Add, expr_tmp, core_expr);
		Stmt my_true_stmt = Move::make(expr_tmp, eadd, MoveType::MemToMem);

		Expr empty = IntImm::make(Type::int_scalar(32), 0);

		Stmt elif = (idv.size() > 0) ? IfThenElse::make(index_based_cond_gen(idv, 0), my_true_stmt, LoopNest::make({}, {})) : my_true_stmt;

		Stmt loop_nest = LoopNest::make(lpv, {elif});
		core_code.push_back(loop_nest);
		return expr_tmp;
	}
}

static void pass_expr(std::vector<std::string> &str, //suffix
					  std::map<std::string, myIndex> &id,
					  std::map<std::string, myVariables> &var,
					  std::map<std::string, myImmediate> &imm,
					  std::vector<Stmt> &core_code){
	Type data_type; //giao
	Expr tmp_zero_expr;
	switch (global_data_type){
	case 0: // float
		data_type = Type::float_scalar(32);
		tmp_zero_expr = FloatImm::make(data_type, 0);
		break;
	case 1:
		data_type = Type::int_scalar(32);
		tmp_zero_expr = IntImm::make(data_type, 0);
		break;
	case 2:
		data_type = Type::uint_scalar(32);
		tmp_zero_expr = UIntImm::make(data_type, 0);
		break;
	default:
		break;
	}
	std::string op_string = "+ - =";
	std::string ase_string = "* / %";
	std::vector<std::string> key; // var_1 var_2 + =

	if (!key.empty())
		key.clear();
	int str_size = str.size();
	for (int i = 0; i < str_size; ++i){
		// std::cout << str[i] << std::endl;
		if (str[i].find("var_") != std::string::npos){
			
			key.push_back(str[i]);
		}
		else if (str[i].find("imm_") != std::string::npos){
			
			key.push_back(str[i]);
		}
		else if (op_string.find(str[i]) != std::string::npos){
			int op_num = 2;
			int pos = key.size() - 1;
			int split = pos;
			bool flag = false;
			int start = pos;
			for (int i = pos; i >= 0; --i){
				if (ase_string.find(key[i]) != std::string::npos){ // * / %
					op_num++;
				}
				else{ // bianliang
					op_num--;
					if (op_num == 1 && flag == false){
						split = i;
						flag = true;
					}
					else if (op_num == 0){
						start = i;
						break;
					}
				}
			}
			std::vector<std::string> one;
			std::vector<std::string> two;
			one.assign(key.begin() + start, key.begin() + split); // the first
			two.assign(key.begin() + split, key.end());

			key.erase(key.begin() + start, key.end());

			bool grad_flag1 = false;
			bool grad_flag2 = false;
			Expr res2 = (str[i][0] == '=') ? var[str[0]].grad_expr : sub_reduction(one, var, id, imm, core_code, data_type, tmp_zero_expr, &grad_flag2);
			Expr res1 = sub_reduction(two, var, id, imm, core_code, data_type, tmp_zero_expr, &grad_flag2);

			myVariables myvar;

			myvar.rw = 0;
			myvar.dim = 1;
			myvar.shape.clear();
			myvar.shape.push_back(0);

			if (str[i][0] == '+'){
				Expr all_res = Binary::make(data_type, BinaryOpType::Add, res2, res1); //
				red_cnt++;
				std::string tmp_red_name = "red_";
				tmp_red_name += std::to_string(red_cnt);
				myvar.name = tmp_red_name;
				myvar.dim = 0;
				myvar.expr = all_res;
				myvar.need_grad = false;
				myvar.grad_expr = tmp_zero_expr;
				myvar.contain_grad = grad_flag1 | grad_flag2;
				var[tmp_red_name] = myvar;
				key.push_back(tmp_red_name);
			}
			else if (str[i][0] == '-'){
				Expr all_res = Binary::make(data_type, BinaryOpType::Sub, res2, res1); //
				red_cnt++;
				std::string tmp_red_name = "red_";
				tmp_red_name += std::to_string(red_cnt);
				myvar.name = tmp_red_name;
				myvar.expr = all_res;
				myvar.dim = 0;
				myvar.need_grad = false;
				myvar.grad_expr = tmp_zero_expr;	
				myvar.contain_grad = grad_flag1 | grad_flag2;
				var[tmp_red_name] = myvar;
				key.push_back(tmp_red_name);
			}
			else if (str[i][0] == '='){

				Stmt all_stmt_2 = Move::make(res2, res1, MoveType::MemToMem);
				core_code.push_back(all_stmt_2);
			}
		}
		else{
			key.push_back(str[i]);
		}
		//A[i,j] = B[i,k] * D[k,k] * E[i,p] * C[p,j] + D[i,j] * E[i,j] * alpha
		//A [B D * C *] [D E * ALPHA *] + =
		//  [B D * C *]
		// A tmp_ =
		//
	}

	return;
}

static Stmt build_tree(std::vector<std::string> &str, //suffix
					   std::map<std::string, myIndex> &id,
					   std::map<std::string, myVariables> &var,
					   std::map<std::string, myImmediate> &imm){
	// const definition

	// index domain, use both_symbol
	std::vector<Stmt> core_code;
	core_code.clear();
	// std::cout << "create finished!" <<std::endl;
	pass_expr(str, id, var, imm, core_code);

	std::vector<Expr> da_loop;
	// for(auto & it : id){
	// 	std::cout << it.second.name <<" ";
	// }
	std::cout << std::endl;
	std::set<std::string> lsset;
	lsset.clear();
	// for (auto &iit : left_symbol){
	// 	for (auto &it : id){
	// 		if (it.second.name.compare(iit) == 0 && lsset.find(iit) == lsset.end() && it.second.ub ){
	// 			// std::cout << it.second.name << std::endl;
	// 			da_loop.push_back(it.second.slf);
	// 			lsset.insert(iit);
	// 		}
	// 	}
	// }
	for(auto &it : var[str[0]].indexes){
		da_loop.push_back(id[it].slf);
	}

	std::vector<myIndex> allif;
	allif.clear();

	for (auto &it : id){
		bool flag = true;

		if (it.second.symbols.size() == 0)
			continue;
		if (it.second.for_loop > 0)
			continue;
		for (auto &iit : it.second.symbols){
			if (left_symbol.find(iit) == left_symbol.end()){
				flag = false;
				break;
			}
		}
		if (flag){
			allif.push_back(it.second);
		}
	}

	Stmt elif = (allif.size() > 0) ? IfThenElse::make(index_based_cond_gen(allif, 0), LoopNest::make({}, core_code), LoopNest::make({}, {})) : LoopNest::make({}, core_code);

	Stmt main_stmt = LoopNest::make(da_loop, {elif});

	return main_stmt;
}

static void update_index_symbols(std::map<std::string, myIndex> &id){
	for (auto &it : id){
		it.second.symbols.clear();
		std::string name = it.second.name;
		std::regex pattern("[a-zA-Z_][a-zA-Z0-9_]*");
		std::smatch result;
		std::string::const_iterator iterStart = name.begin(); //i+j+2
		std::string::const_iterator iterEnd = name.end();
		std::string temp;
		while (std::regex_search(iterStart, iterEnd, result, pattern)){
			temp = result[0];
			it.second.symbols.insert(temp);
			iterStart = result[0].second;
		}
	}
}

static void create_io(std::map<std::string, para> all_name,
					  std::map<std::string, myVariables> &var, std::vector<Expr> &inputs, 
					  std::vector<Expr> &outputs, Json::Value json_object){

	Type index_type = Type::int_scalar(32);
	Type data_type; //giao
	switch (global_data_type){
	case 0: // float
		data_type = Type::float_scalar(32);
		break;
	case 1:
		data_type = Type::int_scalar(32);
		break;
	case 2:
		data_type = Type::uint_scalar(32);
		break;
	default:
		break;
	}

	inputs.clear();
	outputs.clear();

	for(int k = 0; k < json_object["ins"].size(); ++k){
		std::string target_name = json_object["ins"][k].asString();
		std::vector<Expr> idv;
		std::vector<long unsigned int> myshape;
		idv.clear();
		myshape.clear();
		std::string name;
		if(all_name.find(target_name) == all_name.end()) continue;
		for (int i = 0; i < all_name[target_name].shape.size(); ++i){
			int cur_shape = all_name[target_name].shape[i];
			Expr tmp_dom = Dom::make(index_type, 0, (int)all_name[target_name].shape[i]);
			Expr tmp_index = Index::make(index_type, "ti", tmp_dom, IndexType::Spatial);
			idv.push_back(tmp_index);
			myshape.push_back(all_name[target_name].shape[i]);
		}
		Expr var_set = Var::make(data_type, target_name, idv, myshape);
		inputs.push_back(var_set);
	}

	for(int k = 0; k < json_object["outs"].size(); ++k){
		std::string target_name = "d" + json_object["outs"][k].asString();
		std::vector<Expr> idv;
		std::vector<long unsigned int> myshape;
		idv.clear();
		myshape.clear();
		std::string name;
		if(all_name.find(target_name) == all_name.end()) continue;
		for (int i = 0; i < all_name[target_name].shape.size(); ++i){
			Expr tmp_dom = Dom::make(index_type, 0, (int)all_name[target_name].shape[i]);
			Expr tmp_index = Index::make(index_type, "ti", tmp_dom, IndexType::Spatial);
			idv.push_back(tmp_index);
			myshape.push_back(all_name[target_name].shape[i]);
		}
		Expr var_set = Var::make(data_type, target_name, idv, myshape);
		inputs.push_back(var_set);
	}
	for(int k = 0; k < json_object["grad_to"].size(); ++k){
		std::string target_name = "d" + json_object["grad_to"][k].asString();
		std::vector<Expr> idv;
		std::vector<long unsigned int> myshape;
		idv.clear();
		myshape.clear();
		std::string name;
		if(all_name.find(target_name) == all_name.end()) continue;
		for (int i = 0; i < all_name[target_name].shape.size(); ++i){
			Expr tmp_dom = Dom::make(index_type, 0, (int)all_name[target_name].shape[i]);
			Expr tmp_index = Index::make(index_type, "ti", tmp_dom, IndexType::Spatial);
			idv.push_back(tmp_index);
			myshape.push_back(all_name[target_name].shape[i]);
		}
		Expr var_set = Var::make(data_type, target_name, idv, myshape);
		inputs.push_back(var_set);
	}
	// std::cout << "end create io" << std::endl;

}

static void sub_extend(std::vector<std::string>& one, std::map<std::string, myVariables> &var,
					std::map<std::string, std::vector<std::string> > done_map){
	std::vector<std::string> local;

	if(one[0].find("imm_") != std::string::npos){
		return;
	}
	if(one[0].find("DONE") != std::string::npos){
		one.insert(one.begin()+1, done_map[one[0]].begin(), done_map[one[0]].end());
		one.erase(one.begin(), one.begin()+1);
		return;
	}
	local.clear();
	int local_step = 0;
	int length = one.size();

	for(int i = 0; i < length; ++i){
		std::vector<std::string> tmp;
		tmp.clear();
		tmp.insert(tmp.end(), one.begin(), one.end());
		if(one[i][0] == 'v'){
			// generate a new variable
			// var_2 var_3 * var_4 * => var_2 var_6 * var_4 * +
			if(grad_set.find(var[one[i]].name) != grad_set.end()){
				local_step += 1;
				myVariables tmp_var;
				tmp_var.Copy(var[one[i]]);
				tmp_var.need_grad = true;
				var_cnt++;												// expr cnt num
				std::string var_key = "var_" + std::to_string(var_cnt); //A[i] expr1[]
				var[var_key] = tmp_var;
				tmp[i] = var_key;			
				local.insert(local.end(), tmp.begin(), tmp.end());
				if(i != 0){
					local.push_back("+\0");
				}
			}

		}
	}
// update max_step;
	one.clear();
	one.assign(local.begin(), local.end());
	max_step = std::max(max_step, local_step);
}


void ir_extend(std::vector<std::string>& str, std::map<std::string, myVariables> &var,
			   std::map<std::string, myIndex> &id,
			   std::map<std::string, myImmediate> &imm){
	std::vector<std::string> out_str;
	std::vector<bool> flag;
	out_str.clear();
	std::string op_string = "+ - =";
	std::string ase_string = "* / %";
	std::vector<std::string> key; // var_1 var_2 + =

	int done_cnt = 0;
	std::map<std::string, std::vector<std::string> > done_map;

	// var_1 var_2 var_3 + var_4 + imm_1 / =
	// var_1 var_2 
	if (!key.empty())	key.clear();

	int str_size = str.size();
	for (int i = 0; i < str_size; ++i){
		if (str[i].find("var_") != std::string::npos){
			key.push_back(str[i]);
		}
		else if (str[i].find("imm_") != std::string::npos){
			key.push_back(str[i]);
		}
		else if (op_string.find(str[i]) != std::string::npos){
			int op_num = 2;
			int pos = key.size() - 1;
			int split = pos;
			bool flag = false;
			int start = pos;
			for (int i = pos; i >= 0; --i){
				if (ase_string.find(key[i]) != std::string::npos){ // * / %
					op_num++;
				}
				else{ // bianliang
					op_num--;
					if (op_num == 1 && flag == false){
						split = i;
						flag = true;
					}
					else if (op_num == 0){
						start = i;
						break;
					}
				}
			}
			std::vector<std::string> one;
			std::vector<std::string> two;
			one.assign(key.begin() + start, key.begin() + split); // the first
			two.assign(key.begin() + split, key.end());

			key.erase(key.begin() + start, key.end());
			// std::cout <<" one" << std::endl;
			// for(auto & it : one){
			// 	std::cout << it << ' ';
			// }
			// std::cout << std::endl;
			// std::cout <<" two" << std::endl;
			// for(auto & it : two){
			// 	std::cout << it << ' ';
			// }
			// std::cout << std::endl;


			sub_extend(one, var, done_map);
			sub_extend(two, var, done_map);

			if (str[i][0] == '+'){
				done_cnt ++;
				std::string str_cnt = "DONE" + std::to_string(done_cnt);
				done_map[str_cnt].assign(one.begin(),one.end());
				done_map[str_cnt].insert(done_map[str_cnt].end(),two.begin(),two.end());
				done_map[str_cnt].push_back("+\0");
				key.push_back(str_cnt);
			}
			else if (str[i][0] == '-'){
				done_cnt ++;
				std::string str_cnt = "DONE" + std::to_string(done_cnt);
				done_map[str_cnt].assign(one.begin(),one.end());
				done_map[str_cnt].insert(done_map[str_cnt].end(),two.begin(),two.end());
				done_map[str_cnt].push_back("-\0");
				key.push_back(str_cnt);
			}
			else if (str[i][0] == '='){
				out_str.insert(out_str.end(),one.begin(),one.end());	// val  =
				out_str.insert(out_str.end(),two.begin(),two.end());	// DONE
				out_str.push_back("=\0");
				// key.push_back("DONE");
			}
		}
		else{
			key.push_back(str[i]);
		}

	}
	for(auto& it: out_str){
		if(it[0] == 'v'){
			if(var[it].need_grad == true)	std::cout << "d";
			std::cout << var[it].name << " ";}
		else if(it[0] == 'i')
			std::cout << imm[it].name << " ";
		else{
			std::cout << it <<' ';
		}
	}
	std::cout << std::endl;
	for(auto& it: out_str){
		std::cout << it <<' ';
	}
	std::cout << std::endl;
	for(auto& it: out_str){
		if(it[0] != 'v')	continue;
		var[it].Print();
	}
	str.assign(out_str.begin(), out_str.end());
}

std::vector<std::string> match_return_id(std::string reg, std::string target){
	std::regex pattern(reg);
	std::vector<std::string> vec;
	vec.clear();	
	std::smatch result;	
	std::string::const_iterator iterStart = target.begin(); //i+j+2
	std::string::const_iterator iterEnd = target.end();
	std::string temp;
	int i;
	while (std::regex_search(iterStart, iterEnd, result, pattern)){
		temp = result[0];
		vec.push_back(result[0]);
		i = result[0].second - target.begin();
		iterStart = result[0].second;
	}
	return vec;
}
static void broadcast(std::map<std::string, myIndex> &id, std::string old_str, std::string new_str){
	for(auto& it : id){
		std::string s = it.second.name;
		int pos = s.find(old_str);
		// std::cout <<"old_str: " << old_str <<std::endl; 
		if(s.find(old_str) != std::string::npos){
			// std::cout << s << " "<< new_str <<std::endl;
			s.replace(pos, old_str.size(), new_str);
			// std::cout << s << std::endl;
		}
		it.second.name = s;
		it.second.for_loop = -1;
	}
}

struct mod_type{
	std::string s1;	// x/ 15
	std::string s2; // x% 15
	int mod;

};

struct max_type{
	std::string s1;
	int _max;
};

void index_transform(std::vector<std::string>& str, std::map<std::string, myVariables> &var,
			   std::map<std::string, myIndex> &id,
			   std::map<std::string, myImmediate> &imm){
	// pure_access
	bool flag = true;
	std::vector<int> grad_right;
	grad_right.clear();


	for(int i = 1; i < str.size(); ++i){
		if(str[i][0] == 'v' && var[str[i]].need_grad == true){
			grad_right.push_back(i);
			for(auto& it : var[str[i]].indexes){
				if(id[it].for_loop == 0){
					flag = false;
				}
			}
		}
	}
	// std::cout << flag << std::endl;
	if(flag){
		// myVariables tmp_left(var[str[0]]);
		// myVariables tmp_right(var[str[grad_right[0]]]);
		myVariables tmp_left;
		tmp_left.Copy(var[str[0]]);
		myVariables tmp_right;
		tmp_right.Copy(var[str[grad_right[0]]]);

		var[str[0]].Copy(tmp_right);
		for(int i = 0; i < grad_right.size(); ++i){
			var[str[grad_right[i]]].Copy(tmp_left);
		}

		left_symbol.clear();
		left_symbol.insert(var[str[0]].symbols.begin(), var[str[0]].symbols.end());
		var[str[0]].delta_symbols.clear();
		var[str[0]].is_left = true;
		for(int i = 1; i < str.size(); ++i){
			if(str[i][0] != 'v')	continue;
			auto ptr = &var[str[i]];
			ptr->is_left = false;
			ptr->delta_symbols.clear();
			for(auto& itt : ptr->symbols){
				if(left_symbol.find(itt) == left_symbol.end()){
					ptr->delta_symbols.insert(itt);
				}
			}
		}
		
		// std::cout << "after transform" << std::endl;
		// for(auto & it : str){
		// 	if(it[0] != 'v')	continue;
		// 	var[it].Print();
		// }
		// std::cout << "left symbol"<<std::endl;
		// for(auto& it : left_symbol){
		// 	std::cout << it <<' ';
		// }
		// std::cout << std::endl;
	}else{
		//1 regex ("[0-9]+\+|\ * \- -[a-zA-z]+")
		//2 regex ("[a-zA-z]+\+|\-[0-9]+")
		//3 regex ("[a-zA-A]+ / [0-9]")
		// for(int i = 1; i < str.size(); ++i){	
		// 	if(str[i][0] == 'v' && var[str[i]].need_grad == true){
		// 	// grad_right.push_back(i);
		// 	std::cout << var[str[i]].name << ' ';
		// 	for(auto& it : var[str[i]].indexes){
		// 		if(id[it].for_loop == 0){
		// 			std::cout << id[it].name << std::endl;
		// 		}
		// 	}
		// 	}
		// }
		std::map<std::string, max_type> max_map;
		std::map<std::string, mod_type> mod_map;
		max_map.clear();
		mod_map.clear();
		for(int i = 1; i < str.size(); ++i){
			if(str[i][0] != 'v') continue;
			for(int j = 0; j < var[str[i]].dim; ++j){
				if(id[var[str[i]].indexes[j]].for_loop > 0)	continue;
				std::vector<std::string> vec = match_return_id("[a-zA-Z]+", id[var[str[i]].indexes[j]].name);
				// for(auto& itt: vec){
				// 	std::cout << itt <<' ';
				// }
				int state = 0;
				if(id[var[str[i]].indexes[j]].name.find("+") != std::string::npos){
						state = 1;
				}else if(id[var[str[i]].indexes[j]].name.find("*") != std::string::npos){
					state = 2;
				}else if(id[var[str[i]].indexes[j]].name.find("/") != std::string::npos){
					state = 3;
				}else if(id[var[str[i]].indexes[j]].name.find("%") != std::string::npos){
					state = 4;
				}
				if(vec.size() == 2){
					if(state == 1){
						cid_cnt += 1;
						myIndex tmp_id;
						std::string cid_key = "z_" + std::to_string(cid_cnt);
						tmp_id.name = cid_key;
						tmp_id.ub = std::min(tmp_id.ub, var[str[i]].shape[j]);
						tmp_id.lb = std::max(tmp_id.lb, (unsigned long)0);
						tmp_id.for_loop = 1;
						idx_cnt++; // expr cnt num
						std::string idx_key = "idx_" + std::to_string(idx_cnt);
						id[idx_key] = tmp_id; // add to the index map
						var[str[i]].indexes[j] = idx_key;
						std::string old_str = vec[0];
						std::string new_str = "(" + cid_key + "-" + vec[1] + ")";
						broadcast(id, old_str, new_str);
					}
				}else if(vec.size() == 1){
					
					std::vector<std::string> vec2 = match_return_id("[0-9]+", id[var[str[i]].indexes[j]].name);

					
					
					if(state == 1){
						if(max_map.find(vec[0]) == max_map.end())	max_map[vec[0]] = {"", -1};
						int integer = 0;
						if(vec2.size() == 0 || max_map[vec[0]]._max < atoi(vec2[0].c_str())){
							if(vec2.size() == 1)	integer = atoi(vec2[0].c_str());
							cid_cnt += 1;
							myIndex tmp_id;
							std::string cid_key = "t_" + std::to_string(cid_cnt);
							tmp_id.name = cid_key;
							tmp_id.ub = std::min(tmp_id.ub, var[str[i]].shape[j]);
							tmp_id.lb = std::max(tmp_id.lb, (unsigned long)0);
							tmp_id.for_loop = 1;
							idx_cnt++; // expr cnt num
							std::string idx_key = "idx_" + std::to_string(idx_cnt);
							id[idx_key] = tmp_id; // add to the index map
							max_map[vec[0]].s1 = idx_key;
						}
						max_map[vec[0]]._max = std::max(max_map[vec[0]]._max, integer);
					}
					
					
					if(state == 3 || state == 4){
						if(mod_map.find(vec[0]+"_"+vec2[0]) == mod_map.end()){
							mod_map[vec[0]+"_"+vec2[0]] = {"","",atoi(vec2[0].c_str())};
						}
						if(state == 3){
							if(mod_map[vec[0]+"_"+vec2[0]].s1.size() != 0){
								var[str[i]].indexes[j] = mod_map[vec[0]+"_"+vec2[0]].s1;
							}else{
								cid_cnt ++;
								myIndex tmp_id;
								std::string cid_key = "z_" + std::to_string(cid_cnt);
								tmp_id.name = cid_key;
								tmp_id.ub = std::min(tmp_id.ub, var[str[i]].shape[j]);
								tmp_id.lb = std::max(tmp_id.lb, (unsigned long)0);
								tmp_id.for_loop = 1;
								idx_cnt++; // expr cnt num
								std::string idx_key = "idx_" + std::to_string(idx_cnt);
								id[idx_key] = tmp_id; // add to the index map
								var[str[i]].indexes[j] = idx_key;
								mod_map[vec[0]+"_"+vec2[0]].s1 = idx_key;
							}
						}else{
							if(mod_map[vec[0]+"_"+vec2[0]].s2.size() != 0){
								var[str[i]].indexes[j] = mod_map[vec[0]+"_"+vec2[0]].s2;
							}else{
								cid_cnt ++;
								myIndex tmp_id;
								std::string cid_key = "z_" + std::to_string(cid_cnt);
								tmp_id.name = cid_key;
								tmp_id.ub = std::min(tmp_id.ub, var[str[i]].shape[j]);
								tmp_id.lb = std::max(tmp_id.lb, (unsigned long)0);
								tmp_id.for_loop = 1;
								idx_cnt++; // expr cnt num
								std::string idx_key = "idx_" + std::to_string(idx_cnt);
								id[idx_key] = tmp_id; // add to the index map
								var[str[i]].indexes[j] = idx_key;
								mod_map[vec[0]+"_"+vec2[0]].s2 = idx_key;
							}
						}
						
					}
				}
			}
		}
	// mod_map related update
		for(auto & itt : mod_map){
			std::string tmp_s = itt.first;
			std::cout <<"zsb2 " << tmp_s.find("_") << std::endl;
			std::string tmp_old_s1 = tmp_s;
			
			tmp_old_s1[tmp_s.find("_")] = '/';
			std::cout <<"zsb2 " << tmp_old_s1 << std::endl;
			std::string tmp_old_s2 = tmp_s;
			tmp_old_s2[tmp_s.find("_")] = '%';
			std::cout <<"zsb " << tmp_old_s2 << std::endl;
			broadcast(id, tmp_old_s1, id[itt.second.s1].name);
			broadcast(id, tmp_old_s2, id[itt.second.s2].name);
			broadcast(id, tmp_s.substr(0, tmp_s.find("_")), "("+id[itt.second.s1].name + "*" + std::to_string(itt.second.mod) + "+" + id[itt.second.s2].name +")");
		}
	// max_map related update
		std::string max_left = str[grad_right[0]];
		for(auto & it : max_map){
			std::string ost = it.first + "+" + std::to_string(it.second._max);
			// std::cout << "ost_: "<<ost << std::endl;
			for(auto & itt : str){
				if(itt[0] != 'v') continue;
				for(int i = 0; i < var[itt].dim; ++i){
					if(id[var[itt].indexes[i]].name.compare(ost) == 0){
						// std::cout << "find one! " << var[itt].name << std::endl;
						max_left = itt;
						// std::cout << id[it.second.s1].name << std::endl;
						var[itt].indexes[i] = it.second.s1;
					}
				}
			}
		}
	//	max_map broadcast
		bool multiple_access = false;
		for(auto & it : max_map){
			multiple_access = true;
			broadcast(id, it.first, "("+ id[it.second.s1].name + "-"+std::to_string(it.second._max) + ")");
		}
		for(auto& it : str){
			if(it[0] != 'v')	continue;
			var[it].symbols.clear();
			var[it].update_symbols(id);
		}


		if(multiple_access == false){
			myVariables tmp_left;
			tmp_left.Copy(var[str[0]]);
			myVariables tmp_right;
			tmp_right.Copy(var[max_left]);

			var[str[0]].Copy(tmp_right);
			for(int i = 0; i < grad_right.size(); ++i){
				var[str[grad_right[i]]].Copy(tmp_left);
			}
		}else{
			myVariables tmp_left;
			tmp_left.Copy(var[str[0]]);
			myVariables tmp_right;
			tmp_right.Copy(var[max_left]);

			std::vector<std::string> tmp_indexes1 = var[str[0]].indexes;
			var[str[0]].Copy(tmp_right);
			for(int i = 0; i < var[str[0]].dim; ++i){
				auto ptr = &id[var[str[0]].indexes[i]];
				myIndex tmpId;
				tmpId.name = ptr->name;
				tmpId.lb = 0;
				tmpId.ub = var[str[grad_right[0]]].shape[i];
				tmpId.for_loop = 1;
				idx_cnt++; // expr cnt num
				std::string idx_key = "idx_" + std::to_string(idx_cnt);
				id[idx_key] = tmpId; // add to the index map
				tmp_indexes1[i] = idx_key;
			}
			
			var[str[0]].indexes.assign(tmp_indexes1.begin(), tmp_indexes1.end());



			for(int i = 0; i < grad_right.size(); ++i){
				std::vector<std::string> tmp_indexes = var[str[grad_right[i]]].indexes;
				var[str[grad_right[i]]].Copy(tmp_left);
				for(int j = 0; j < tmp_indexes.size(); ++j){
					id[tmp_indexes[j]].ub = id[var[str[grad_right[i]]].indexes[j]].ub;
				}

				var[str[grad_right[i]]].indexes.assign(tmp_indexes.begin(), tmp_indexes.end());
			}
		}


		for(auto& it : str){
			if(it[0] != 'v')	continue;
			var[it].symbols.clear();
			var[it].update_symbols(id);
		}

		left_symbol.clear();
		left_symbol.insert(var[str[0]].symbols.begin(), var[str[0]].symbols.end());
		var[str[0]].delta_symbols.clear();
		var[str[0]].is_left = true;
		for(int i = 1; i < str.size(); ++i){
			if(str[i][0] != 'v')	continue;
			auto ptr = &var[str[i]];
			ptr->is_left = false;
			ptr->delta_symbols.clear();
			for(auto& itt : ptr->symbols){
				if(left_symbol.find(itt) == left_symbol.end()){
					ptr->delta_symbols.insert(itt);
				}
			}
		}
		
	}
			std::cout << "after transform" << std::endl;

		
		for(auto & it : str){
			if(it[0] != 'v')	std::cout << it;
			// var[it].Print();
			else{
				std::cout << var[it].name;
			for(int i = 0; i < var[it].dim; ++i){
				std::cout << "[" <<id[var[it].indexes[i]].lb << ' '<< id[var[it].indexes[i]].name << " " << id[var[it].indexes[i]].ub << "]";
			}
			}
			std::cout << " ";
		}
		std::cout <<std::endl;

				for(auto & it : str){
			if(it[0] != 'v')	std::cout << it;
			// var[it].Print();
			else{
				std::cout << var[it].name;
			for(int i = 0; i < var[it].dim; ++i){
				std::cout << "[" << var[it].indexes[i] << "]";
			}
			}
			std::cout << " ";
		}
		std::cout <<std::endl;
		// std::cout << "left symbol"<<std::endl;
		// for(auto& it : left_symbol){
		// 	std::cout << it <<' ';
		// }
		std::cout << std::endl;

}





void process(std::string sin, std::string sout){

	Json::Reader myreader;
	Json::Value json_object;

	std::map<std::string, myVariables> var; //map['A'] = {}
	std::map<std::string, myIndex> index;	//map["i+j+2"] = {}
	std::map<std::string, myImmediate> imm; //map["tmp1"] = {}

	var.clear();
	index.clear();
	imm.clear();
	grad_set.clear();
	std::ifstream ifs(sin, std::ios::binary);

	// json parser
	if (!myreader.parse(ifs, json_object)){ // A[] = B+C;
		//ifs.close();
		std::cout << "open error" << std::endl;
		return;
	}
	ifs.close();
	int grad_num = json_object["grad_to"].size();
	std::vector<Stmt> all_main_stmt;
	all_main_stmt.clear();
	std::map<std::string, para> all_name;
	all_name.clear();
	for(int i = 0; i < grad_num; ++i){
		var.clear();
		index.clear();
		imm.clear();
		grad_set.clear();
		process_io(json_object, var, grad_set, i);
		std::vector<std::string> kernels = split_kernels(json_object); // A = B+C; C = B+A;

		max_step = 1;
		if (!left_symbol.empty())
			left_symbol.clear();
		if (!both_symbol.empty())
			both_symbol.clear();
		std::string myir = process_kernel(kernels[0], var, index, imm);
		update_index_symbols(index);

		std::vector<std::string> suffix = inver(myir);
		
		ir_extend(suffix, var, index, imm);
		index_transform(suffix, var, index, imm);
		
		for(auto & it : suffix){
			if(it[0] != 'v')	continue;
			if(var[it].need_grad)	all_name["d"+var[it].name] = {var[it].need_grad, var[it].shape};
			else all_name[var[it].name] = {var[it].need_grad, var[it].shape};
		}
		create_expr(var, imm, index);
		all_main_stmt.push_back(build_tree(suffix, index, var, imm));

		left_symbol.clear();
		both_symbol.clear();
		index.clear();
		imm.clear();
	}
	std::vector<Expr> all_inputs;
	std::vector<Expr> all_outputs;
	create_io(all_name, var, all_inputs, all_outputs, json_object);
	Group kernel = Kernel::make(json_object["name"].asString(), all_inputs, {}, all_main_stmt, KernelType::CPU);

	IRVisitor visitor;
	kernel.visit_group(&visitor);
	// mutator
	IRMutator mutator;
	kernel = mutator.mutate(kernel);

	// printer
	IRPrinter printer;
	std::time_t result = std::time(NULL);
	std::string code = "#include \"../run2.h\" \n" + printer.print(kernel);
	std::ofstream ofile(sout, std::ios::out);
	ofile << "//Produced at : ";
	ofile << std::asctime(std::localtime(&result));
	ofile << code;
	ofile << "//Finished at : ";
	ofile << std::asctime(std::localtime(&result));
	ofile.close();
}

int main(){
	std::string in_prefix = "./cases/case";
	std::string out_prefix = "./kernels/grad_case";
	for (int i = 1; i <= 10; ++i){
		var_cnt = 0;
		idx_cnt = 0;
		imm_cnt = 0;
		tmp_cnt = 0;
		red_cnt = 0;
		cid_cnt = 0;
		process(in_prefix + std::to_string(i) + ".json", out_prefix + std::to_string(i) + ".cc");
	}
}
