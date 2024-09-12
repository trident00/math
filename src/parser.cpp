#include <iostream>
#include <stdio.h>
#include <array>
#include <utility>
#include <cmath>
#include <variant>
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <set>
#include <stack>
#include <queue>
#include <cassert>
#include "parser.h"
#include "ast.h"


using String = std::string;

#define u64		unsigned long long int

#define PI		3.141592653589793238462643383279502884197169399375105820974944592307816	

std::unordered_map<String, TokenType> token_map = {
	{"(", TokenType::LPAREN},
	{")", TokenType::RPAREN},
	{"+", TokenType::OPERATOR_ADD},
	{"-", TokenType::OPERATOR_SUBTRACT},
	{"*", TokenType::OPERATOR_MULTIPLY},
	{"/", TokenType::OPERATOR_DIVIDE},
	{"^", TokenType::OPERATOR_POWER},

	{"sin", TokenType::FUNCTION_SIN},
	{"cos", TokenType::FUNCTION_COS},
	{"tan", TokenType::FUNCTION_TAN},
	{"csc", TokenType::FUNCTION_CSC},
	{"sec", TokenType::FUNCTION_SEC},
	{"cot", TokenType::FUNCTION_COT},
	{"abs", TokenType::FUNCTION_ABS},
	{"sqrt", TokenType::FUNCTION_SQRT},
	{"log", TokenType::FUNCTION_LOG},

	{"e", TokenType::NUMBER_E},
	{"x", TokenType::VARIABLE_X},
	{"y", TokenType::VARIABLE_Y},
	{"z", TokenType::VARIABLE_Z},
	{"pi", TokenType::NUMBER_PI}
};

bool is_function(Token* t) { return (t->type >= FUNCTION_SIN && t->type <= FUNCTION_LOG); }
bool is_variable(Token* t) { return (t->type >= VARIABLE_X && t->type <= VARIABLE_Z); }
bool is_operator(Token* t) { return (t->type >= OPERATOR_ADD && t->type <= OPERATOR_POWER); }
bool is_paren(Token* t) { return t->type == LPAREN || t->type == RPAREN; }
bool is_number(Token* t) { return (t->type >= NUMBER && t->type <= NUMBER_E); }

Token* create_token_mult()
{
	Token* mult = new Token();
	mult->value = "*";
	mult->type = OPERATOR_MULTIPLY;
	return mult;
}

Token* create_token_add()
{
	Token* t = new Token();
	t->value = "+";
	t->type = OPERATOR_ADD;
	return t;
}

Token* create_token_negative()
{
	Token* negative = new Token();
	negative->num = -1;
	negative->type = NUMBER;
	return negative;
}

std::vector<Token*> lexer(const String& input)
{
	std::vector<Token*> result;
	String num_buffer;
	String str_buffer;
	for (u64 i = 0; i < input.length(); i++) {
		char c = input[i];
		if (isdigit(c) || c == '.') {
			num_buffer += c;
			if (!isdigit(input[i+1]) && input[i+1] != '.') {
				Token* n_token = new Token();
				n_token->type = NUMBER;
				// Handle negative numbers:
				if (!result.empty() && result.back()->type == OPERATOR_SUBTRACT) {
					result.pop_back();
					n_token->num = -std::stod(num_buffer);
					result.emplace_back(create_token_add());
				} else {
					n_token->num = std::stod(num_buffer);
				}
				result.emplace_back(n_token);
				num_buffer.clear();
			}
		} else if (c != ' ') {
			str_buffer += c;
			auto it = token_map.find(str_buffer);
			if (it != token_map.end()) {
				Token* f_token = new Token();
				f_token->value = str_buffer;
				f_token->type = it->second;
				if (f_token->type == NUMBER_PI) f_token->num = PI;
				if (f_token->type == NUMBER_E) f_token->num = E;
				// Handle implicit mult: ie. 3x = 3*x
				if (!result.empty() && !is_operator(f_token)) {
					Token* adj = result.back();
					if (f_token->type == LPAREN && adj->type == RPAREN) {
						// 'Binomial' expression: (x+y)(x+y) = (x+y)*(x+y)
						result.emplace_back(create_token_mult());
					} else if ((is_number(adj) || is_variable(adj))
							&& (is_function(f_token) || f_token->type == LPAREN || is_variable(f_token))) {
						// Number/Variable before Paren/Function/Variable: 3(x) = 3*(x), 3sin(x) = 3*sin(x)
						result.emplace_back(create_token_mult());
					} else if ((is_variable(f_token) || is_function(f_token)) && adj->type == OPERATOR_SUBTRACT) {
						result.pop_back();
						result.emplace_back(create_token_add());
						result.emplace_back(create_token_negative());
						result.emplace_back(create_token_mult());
					}
				}
				result.emplace_back(f_token);
				str_buffer.clear();	
			}
		}
	}
	// @DEBUG
	For (input) {
		std::cout<<item<<" ";
	}
	std::cout<<std::endl;
	For (result) {
		std::cout<<print_token(item)<<" ";
	}
	std::cout<<std::endl;
	return result;
}

struct Operator
{
	TokenType type;
	int precedence;
	bool associativity; // true is left
};

std::unordered_map<TokenType, Operator> operator_map = {
	{TokenType::OPERATOR_POWER, {TokenType::OPERATOR_POWER, 4, false}},
	{TokenType::OPERATOR_MULTIPLY, {TokenType::OPERATOR_MULTIPLY, 3, true}},
	{TokenType::OPERATOR_DIVIDE, {TokenType::OPERATOR_DIVIDE, 3, true}},
	{TokenType::OPERATOR_ADD, {TokenType::OPERATOR_ADD, 2, true}},
	{TokenType::OPERATOR_SUBTRACT, {TokenType::OPERATOR_SUBTRACT, 2, true}}

};

std::queue<Token*> parser(std::vector<Token*> tokens)
{
	std::queue<Token*> output;
	std::stack<Token*> operators;


	// Shunting yard algorithm, see https://en.wikipedia.org/wiki/Shunting_yard_algorithm
	for (Token* token : tokens) {

		// Numbers go straight to output
		if (is_number(token)) {
			output.push(token);
		} else if (is_variable(token)) {
			output.push(token);
		} else if (is_function(token)) {
			operators.push(token);
		} else if (is_operator(token)) {
			auto o1_it = operator_map.find(token->type);
			Operator o_1 = o1_it->second;

			while (!operators.empty() && is_operator(operators.top()) && operators.top()->type != TokenType::LPAREN) {
				auto o2_it = operator_map.find(operators.top()->type);
				Operator o_2 = o2_it->second;
				if (o_2.precedence > o_1.precedence || (o_1.precedence == o_2.precedence && o_1.associativity)) {
					output.push(operators.top());
					operators.pop();
				} else break;
			}
			operators.push(token);
		} else if (token->type == TokenType::LPAREN) {
			operators.push(token);
		} else if (token->type == TokenType::RPAREN) {
			while (!operators.empty() && operators.top()->type != TokenType::LPAREN) {
				assert(!operators.empty()); // missing matching parenthesis
				output.push(operators.top());
				operators.pop();
			}
			assert(operators.top()->type == TokenType::LPAREN); // missing matching parenthesis
			operators.pop();
			if (!operators.empty() && is_function(operators.top())) {
				output.push(operators.top());
				operators.pop();
			}
		}
	}

	while (!operators.empty()) {
		assert(operators.top()->type != TokenType::LPAREN); // missing matching parenthesis
		output.push(operators.top());
		operators.pop();
	}

	print_queue(output);

	return output;
}

std::queue<Token*> tokenize(String f)
{
	return parser(lexer(f));
}



AstNode* to_ast(std::queue<Token*> tokens)
{
	std::stack<AstNode*> buffer_stack;
	while (!tokens.empty()) {
		Token* t = tokens.front();
		std::cout << "Processing Token: " << print_token(t) << std::endl;
		if (is_number(t)) {
			buffer_stack.push(create_node_number(t));
			tokens.pop();
		} else if (is_operator(t)) {
			assert(buffer_stack.size() >= 2);
			AstNode* right = buffer_stack.top();
			buffer_stack.pop();
			AstNode* left = buffer_stack.top();
			buffer_stack.pop();
			buffer_stack.push(create_node_operator(t, left, right));
			tokens.pop();
		} else if (is_variable(t)) {
			buffer_stack.push(create_node_variable(t));
			tokens.pop();
		} else if (is_function(t)) {
			assert(buffer_stack.size() >= 1);
			AstNode* arg = buffer_stack.top();
			buffer_stack.pop();
			buffer_stack.push(create_node_function(t, arg));
			tokens.pop();
		}
		//std::cout << " -> Buffer Stack Size: " << buffer_stack.size() << std::endl;
	}
	assert(buffer_stack.size() == 1);
	print_ast_tree(buffer_stack.top());
	return buffer_stack.top();
}

// MATH FUNCTIONS
// BASIC
double power(double a, double b) { return std::pow(a, b); }
double multiply(double a, double b) { return a * b; }
double divide(double a, double b) { if (b == 0); return a / b; }
double add(double a, double b) { return a + b; }
double subtract(double a, double b) { return a - b; }
// TRIG (RADIANS)
double degree(double a) { return a * PI / 180.0; }
double math_sin(double a) { return std::sin(a); }
double math_cos(double a) { return std::cos(a); }
double math_tan(double a) { return std::tan(a); }
double math_csc(double a) { return 1.0 / std::sin(a); }
double math_sec(double a) { return 1.0 / std::cos(a); }
double math_cot(double a) { return 1.0 / std::tan(a); }
// Others
double math_abs(double a) { return a >=0 ? a : -a; }
double math_sqrt(double a) { return power(a, 0.5); }
double math_log(double a) { if (a <= 0); return std::log(a); }

std::unordered_map<TokenType, double(*)(double, double)> binary_operation = {
	{OPERATOR_POWER,		power},
	{OPERATOR_MULTIPLY,		multiply},
	{OPERATOR_DIVIDE,		divide},
	{OPERATOR_ADD, 			add},
	{OPERATOR_SUBTRACT, 	subtract}
};

std::unordered_map<TokenType, double(*)(double)> unary_operation = {
	{FUNCTION_SIN, 			math_sin},
	{FUNCTION_COS, 			math_cos},
	{FUNCTION_TAN, 			math_tan},
	{FUNCTION_CSC, 			math_csc},
	{FUNCTION_SEC, 			math_sec},
	{FUNCTION_COT, 			math_cot},
	{FUNCTION_ABS,			math_abs},
	{FUNCTION_SQRT,			math_sqrt},
	{FUNCTION_LOG, 			math_log},
};


double ast_solve(AstNode* node, double x, double y, double z)
{
	switch (node->type) {
		case NODE_NUMBER:
			return node->number_value;
		case NODE_OPERATOR:
			if (node->OperatorNode.left && node->OperatorNode.right) {
				double left = ast_solve(node->OperatorNode.left,x,y,z);
				double right = ast_solve(node->OperatorNode.right,x,y,z);
				return binary_operation.at(node->OperatorNode.op)(left, right);
			} else throw std::runtime_error("Missing branch(es) for node NODE_OPERATOR somewhere.");
		case NODE_VARIABLE:
			switch (node->variable) {
				case VARIABLE_X:	return x;
				case VARIABLE_Y:	return y;
				case VARIABLE_Z:	return z;
				default: throw std::runtime_error("Non-variable assigned as NODE_VARIABLE.");
			}
		case NODE_FUNCTION:
			if (node->FunctionNode.argument) {
				double arg = ast_solve(node->FunctionNode.argument,x,y,z);
				return unary_operation.at(node->FunctionNode.function)(arg);
			} else throw std::runtime_error("Missing argument for node NODE_FUNCTION.");
		default:	throw std::runtime_error("Unknown NODE in ast_solve");
	}
	return 0.0;
}


AstNode* ast_tree(String input)
{
	return to_ast(tokenize(input));
}

String ast_to_string(AstNode* node)
{
	String result;
	AstNode* current = node;
	while (current) {
		switch (current->type) {
			case NODE_VARIABLE:
			case NODE_NUMBER:
				result += node->symbol;
				current = nullptr;
				break;
			case NODE_OPERATOR:
				result += ast_to_string(node->OperatorNode.left);
				result += node->symbol;
				result += ast_to_string(node->OperatorNode.right);
				current = nullptr;
				break;
			case NODE_FUNCTION:
				result += node->symbol;
				result += "(" + ast_to_string(node->FunctionNode.argument) + ")";
				current = nullptr;
				break;
			default:	throw std::runtime_error("Unknown NODE in ast_to_string");
		}
	}
	return result;
}

void debug()
{
	String inp = "e^x/(1+e^x)";
//	String inp = "2xy(sin(2x)) + (3*x-4) / (x+2)";
	std::cout<< "\n" << inp << std::endl;
	auto s = lexer(inp);
	std::cout<<"\nlexer_main-:\n";
	For (s) {
		std::cout<<print_token(item);
	}
	std::cout<<"\n\nto_ast-:\n";
	AstNode* tree = to_ast(parser(s));
	std::cout<<"\n"<<ast_solve(tree,2.018)<<std::endl;
}

