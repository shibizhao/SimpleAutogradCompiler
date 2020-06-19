/* 
Compiler Project 01: JSON file based C++ code generator.

Author: Bizhao Shi, Shaobo Zhou, Minggui Teng

Email: {shi_bizhao, telsazhou, minggui_teng}@pku.edu.cn

Date: May/17/2020

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
const int typeNum = 4;

std::set<std::string> left_symbol;
std::set<std::string> both_symbol;
int global_data_type;

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
	std::vector<int> shape;			  //1  A[i+j+3, i+j] a[][]
	std::vector<std::string> indexes; // index : idx_t
	std::set<std::string> symbols;	  //i+j+k i,j,k
	std::set<std::string> delta_symbols;
	int rw;
	Expr expr;
	myVariables() {}
	myVariables(std::string s, int data_type, int flag); // flag is read or write sign
	void Print();
	void update_symbols(std::map<std::string, myIndex> &id);
};

class myIndex{
public:
	std::string name; //I+j+2
	int lb;
	int ub;
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
	// 	std::cout << "shape: \n";
	// 	for(int i = 0; i < (int)shape.size(); ++i){
	// 		std::cout << shape[i] << std::endl;
	// 	}
	// //	std::cout << std::endl;
	// 	for(auto & it : symbols){
	// 		std::cout << it << ' ';
	// 	}
	// 	std::cout << std::endl;std::cout << std::endl;
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
	lb = -1;
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

static void process_io(const Json::Value &json_object, std::map<std::string, myVariables> &var){
	global_data_type = judge_type(json_object["data_type"].asString());
	//out
	int out_length = json_object["outs"].size();
	for (int i = 0; i < out_length; ++i){
		std::string tmp_name = json_object["outs"][i].asString();
		myVariables add_var = myVariables(tmp_name, global_data_type, -1);
		add_var.is_left = -1;
		var[tmp_name] = add_var; // read is 0
	}
	// in
	int in_length = json_object["ins"].size();
	for (int i = 0; i < in_length; ++i){
		std::string tmp_name = json_object["ins"][i].asString();
		myVariables add_var = myVariables(tmp_name, global_data_type, 1);
		add_var.is_left = -1;
		var[tmp_name] = add_var; // read is 0
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
			std::string id_name = tmp_index.substr(pos, end_pos - pos);
			std::string root = CheckIsIndex(index, id_name);
			if (root.compare("NOT_FOUND") == 0){
				myIndex id;
				// check the valid
				id.name = id_name;
				id.ub = std::min(id.ub, tmp_var.shape[i]);
				id.lb = std::max(id.lb, 0);
				id.for_loop = id.CheckIndex(id.name);

				idx_cnt++; // expr cnt num
				std::string idx_key = "idx_" + std::to_string(idx_cnt);
				index[idx_key] = id; // add to the index map
				tmp_var.indexes.push_back(idx_key);
			}
			else{ //already exists
				index[root].ub = std::min(index[root].ub, tmp_var.shape[i]);
				index[root].lb = std::max(index[root].lb, 0);
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
		it.second.dom = Dom::make(index_type, it.second.lb, it.second.ub);
		it.second.slf = Index::make(index_type, it.second.name, it.second.dom, IndexType::Spatial);
	}

	Type data_type; //giao
	switch (global_data_type){
	case 0:
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
	for (auto &it : var){
		if (it.first[0] != 'v')
			continue;
		if (it.second.dim == 1 && it.second.shape[0] == 1){
			it.second.expr = Var::make(data_type, it.second.name, {}, {1});
		}
		else{
			std::vector<Expr> tmp_id;
			std::vector<long unsigned int> my_shape;
			for (int dim_i = 0; dim_i < it.second.dim; ++dim_i){
				tmp_id.push_back(id[it.second.indexes[dim_i]].slf);
				my_shape.push_back((long unsigned int)(it.second.shape[dim_i]));
			}
			it.second.expr = Var::make(data_type, it.second.name, tmp_id, my_shape);
		}
	}
	for (auto &it : imm){
		switch (it.second.type){
		case 0:{ //float
			Type imm_type = Type::float_scalar(32);
			it.second.expr = FloatImm::make(imm_type, atof(it.second.value.c_str()));
			break;
		}
		case 1:{
			Type imm_type = Type::int_scalar(32);
			it.second.expr = IntImm::make(imm_type, atoi(it.second.value.c_str()));
			break;
		}
		case 2:{
			Type imm_type = Type::uint_scalar(32);
			it.second.expr = UIntImm::make(imm_type, (unsigned int)atoi(it.second.value.c_str()));
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
						 Type data_type){
	std::stack<Expr> s;

	std::string op_string = "*/%";
	for (auto &it : str){ //var_1 *
		if (op_string.find(it) != std::string::npos){
			Expr top1 = s.top();
			s.pop();
			Expr top2 = s.top();
			s.pop();
			s.push(expr_gen(top2, top1, op_string.find(it)));
		}
		else{
			if (it.find("var_") != std::string::npos){
				s.push(var[it].expr);
			}
			else if (it.find("red_") != std::string::npos){
				s.push(var[it].expr);
			}
			else if (it.find("imm_") != std::string::npos){
				s.push(imm[it].expr);
			}
		}
	}
	Expr top = s.top();
	s.pop();

	return top;
}

static Expr sub_reduction(std::vector<std::string> &str, //suffix
						  std::map<std::string, myVariables> &var,
						  std::map<std::string, myIndex> &id,
						  std::map<std::string, myImmediate> &imm,
						  std::vector<Stmt> &core_code,
						  Type data_type){
	// red_1
	if (str.size() == 1 && str[0].find("red_") != std::string::npos){
		return var[str[0]].expr;
	}
	else{ // var_1 imm_

		Expr core_expr = pass_sub_str(str, id, var, imm, data_type);

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

		if (loop_se.size() == 0)
			return core_expr;

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
				if (var[str[i]].dim == 1 && var[str[i]].shape[0] == 1)
					continue;
				for (auto &it : var[str[i]].indexes){ //idx_1/2/3
					if (id[it].for_loop > 0)
						continue;
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
	std::string op_string = "+ - =";
	std::string ase_string = "* / %";
	std::stack<Expr> s;			  // var_1 var_2
	std::vector<std::string> key; // var_1 var_2 + =
	while (!s.empty())
		s.pop();
	if (!key.empty())
		key.clear();
	int str_size = str.size();
	for (int i = 0; i < str_size; ++i){
		if (str[i].find("var_") != std::string::npos){
			s.push(var[str[i]].expr);
			key.push_back(str[i]);
		}
		else if (str[i].find("imm_") != std::string::npos){
			s.push(imm[str[i]].expr);
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
					s.pop();
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
			Expr res2 = sub_reduction(one, var, id, imm, core_code, data_type);
			Expr res1 = sub_reduction(two, var, id, imm, core_code, data_type);

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
				myvar.expr = all_res;
				var[tmp_red_name] = myvar;
				s.push(all_res);
				key.push_back(tmp_red_name);
			}
			else if (str[i][0] == '-'){
				Expr all_res = Binary::make(data_type, BinaryOpType::Sub, res2, res1); //
				red_cnt++;
				std::string tmp_red_name = "red_";
				tmp_red_name += std::to_string(red_cnt);
				myvar.name = tmp_red_name;
				myvar.expr = all_res;
				var[tmp_red_name] = myvar;
				s.push(all_res);
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
	create_expr(var, imm, id);
	pass_expr(str, id, var, imm, core_code);

	std::vector<Expr> da_loop;
	for (auto &iit : left_symbol){
		for (auto &it : id){
			if (it.second.name.compare(iit) == 0){
				da_loop.push_back(it.second.slf);
			}
		}
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

static void create_io(std::map<std::string, myVariables> &var,
					  std::vector<Expr> &inputs, std::vector<Expr> &outputs){
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

	for (auto &it : var){
		if (it.second.rw == 0)
			continue;
		std::vector<Expr> idv;
		std::vector<long unsigned int> myshape;
		idv.clear();
		myshape.clear();
		for (int i = 0; i < it.second.dim; ++i){
			Expr tmp_dom = Dom::make(index_type, 0, it.second.shape[i]);
			Expr tmp_index = Index::make(index_type, "ti", tmp_dom, IndexType::Spatial);
			idv.push_back(tmp_index);
			myshape.push_back(it.second.shape[i]);
		}
		Expr var_set = Var::make(data_type, it.second.name, idv, myshape);
		if (it.second.rw == 1){
			inputs.push_back(var_set);
		}
		else{
			outputs.push_back(var_set);
		}
	}
	int opt_size = outputs.size();
	for (int i = 0; i < opt_size; ++i){
		inputs.push_back(outputs[i]);
	}
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
	std::ifstream ifs(sin, std::ios::binary);

	// json parser
	if (!myreader.parse(ifs, json_object)){ // A[] = B+C;
		//ifs.close();
		std::cout << "open error" << std::endl;
		return;
	}
	ifs.close();

	process_io(json_object, var);

	std::vector<std::string> kernels = split_kernels(json_object); // A = B+C; C = B+A;

	std::vector<Stmt> all_main_stmt;
	all_main_stmt.clear();
	int nkernel = kernels.size();
	for (int i = 0; i < nkernel; ++i){ // every loop, every block of for-loop is generated
		if (!left_symbol.empty())
			left_symbol.clear();
		if (!both_symbol.empty())
			both_symbol.clear();
		std::string myir = process_kernel(kernels[i], var, index, imm);
		update_index_symbols(index);
		std::vector<std::string> suffix = inver(myir);
		all_main_stmt.push_back(build_tree(suffix, index, var, imm));

		left_symbol.clear();
		both_symbol.clear();
		index.clear();
		imm.clear();
	}

	std::vector<Expr> all_inputs;
	std::vector<Expr> all_outputs;
	create_io(var, all_inputs, all_outputs);

	Group kernel = Kernel::make(json_object["name"].asString(), all_inputs, {}, all_main_stmt, KernelType::CPU);

	IRVisitor visitor;
	kernel.visit_group(&visitor);
	// mutator
	IRMutator mutator;
	kernel = mutator.mutate(kernel);

	// printer
	IRPrinter printer;
	std::time_t result = std::time(NULL);
   
	std::string code = "\n #include \"../run.h\" \n" + printer.print(kernel);

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
	std::string out_prefix = "./kernels/kernel_case";
	process("./cases/example.json", "./kernels/kernel_example.cc");
	for (int i = 1; i <= 10; ++i){
		var_cnt = 0;
		idx_cnt = 0;
		imm_cnt = 0;
		tmp_cnt = 0;
		red_cnt = 0;
		process(in_prefix + std::to_string(i) + ".json", out_prefix + std::to_string(i) + ".cc");
	}
}
