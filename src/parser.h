#pragma once
#include <string>
#include <queue>
#include <unordered_map>

#define For(container) for (auto& item : container)

#define E		2.7182818284590452353602874713527
using String = std::string;
struct AstNode;

enum TokenType
{
	NUMBER,
	NUMBER_PI,
	NUMBER_E,

	VARIABLE_X,
	VARIABLE_Y,
	VARIABLE_Z,
	
	LPAREN,
	RPAREN,
	OPERATOR_ADD,
	OPERATOR_SUBTRACT,
	OPERATOR_MULTIPLY,
	OPERATOR_DIVIDE,
	OPERATOR_POWER,

	// Unary
	FUNCTION_SIN,
	FUNCTION_COS,
	FUNCTION_TAN,
	FUNCTION_CSC,
	FUNCTION_SEC,
	FUNCTION_COT,
	FUNCTION_ABS,
	FUNCTION_SQRT,
	FUNCTION_LOG
};

struct Token
{
	String value; // String value associated with token
	double num; // Number associated with token
	TokenType type; // Type of token

	Token() {}
	~Token() {}
};

bool is_function(Token* t);
bool is_variable(Token* t);
bool is_operator(Token* t); 
bool is_number(Token* t);
std::vector<Token*> lexer(const String &input);
std::queue<Token*> parser(std::vector<Token*> tokens);

std::queue<Token*> tokenize(String f);
AstNode* ast_tree(String input);
double ast_solve(AstNode* node, double x=0, double y=0, double z=0);
String ast_to_string(AstNode* node);
