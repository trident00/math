#pragma once
#include <string>
#include "parser.h"

using String = std::string;
enum NodeType
{
	NODE_NUMBER,
	NODE_VARIABLE,
	NODE_OPERATOR,
	NODE_FUNCTION
};

struct AstNode
{
	NodeType type;
	String symbol; // For printing
	union {
		// NODE_NUMBER
		double number_value;
		// NODE_VARIABLE
		TokenType variable;
		// NODE_OPERATOR
		struct {
			TokenType op;
			AstNode* left;
			AstNode* right;
		} OperatorNode;
		// NODE_FUNCTION
		struct {
			TokenType function;
			AstNode* argument;
		} FunctionNode;
	};

	AstNode() {}	// Constructor
	~AstNode() {}	// Destructor
};

void print_queue(std::queue<Token*> q);
AstNode* create_node_number(Token* token);
AstNode* create_node_variable(Token* token);
AstNode* create_node_operator(Token* token, AstNode* left, AstNode* right);
AstNode* create_node_function(Token* token, AstNode* arg);

AstNode* ast_tree(String input);

// Debug
String to_string(TokenType t);
String to_string(NodeType t);
String print_token(Token* t);
std::string print_ast(AstNode* node);
void print_ast_tree(const AstNode* node, int depth=0, String prefix="", bool first_flag=true);
