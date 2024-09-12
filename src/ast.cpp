#include <iostream>
#include <vector>
#include <variant>
#include <cmath>
#include <unordered_map>
#include <stack>
#include <queue>
#include <cassert>
#include "ast.h"
#include "parser.h"

#define For(container) for (auto& item : container)

#define PI		3.141592653589793238462643383279502884197169399375105820974944592307816	
#define E		2.7182818284590452353602874713527




//
//	Create Node Functions
//
AstNode* create_node_number(Token* token)
{
	AstNode* node = new AstNode();
	node->type = NODE_NUMBER;
	node->number_value = token->num;
	node->symbol = std::to_string(token->num);
	return node;
}

AstNode* create_node_variable(Token* token)
{
	AstNode* node = new AstNode();
	node->type = NODE_VARIABLE;
	node->variable = token->type;
	node->symbol = token->value;
	return node;
}

AstNode* create_node_operator(Token* token, AstNode* left, AstNode* right)
{
	AstNode* node = new AstNode();
	node->type = NODE_OPERATOR;
	node->OperatorNode.op = token->type;
	node->OperatorNode.left = left;
	node->OperatorNode.right = right;
	node->symbol = token->value;
	return node;
}

AstNode* create_node_function(Token* token, AstNode* arg)
{
	AstNode* node = new AstNode();
	node->type = NODE_FUNCTION;
	node->FunctionNode.function = token->type;
	node->FunctionNode.argument = arg;
	node->symbol = token->value;
	return node;
}

//
//	DEBUG
//
String to_string(TokenType t)
{
    switch (t) {
		case NUMBER:				return "NUMBER";
		case NUMBER_PI:				return "NUMBER_PI";
		case NUMBER_E:				return "NUMBER_E";
		case VARIABLE_X:			return "VARIABLE_X";
		case VARIABLE_Y:			return "VARIABLE_Y";
		case VARIABLE_Z:			return "VARIABLE_Z";
		case LPAREN:				return "LPAREN";
		case RPAREN:				return "RPAREN";
		case OPERATOR_ADD:			return "OPERATOR_ADD";
		case OPERATOR_SUBTRACT:		return "OPERATOR_SUBTRACT";
		case OPERATOR_MULTIPLY:		return "OPERATOR_MULTIPLY";
		case OPERATOR_DIVIDE:		return "OPERATOR_DIVIDE";
		case OPERATOR_POWER:		return "OPERATOR_POWER";
		case FUNCTION_SIN:			return "FUNCTION_SIN";
		case FUNCTION_COS:			return "FUNCTION_COS";
		case FUNCTION_TAN:			return "FUNCTION_TAN";
		case FUNCTION_CSC:			return "FUNCTION_CSC";
		case FUNCTION_SEC:			return "FUNCTION_SEC";
		case FUNCTION_COT:			return "FUNCTION_COT";
		case FUNCTION_LOG:			return "FUNCTION_LOG";
		default:		return "UNKNOWN";
    }
}

String to_string(NodeType t)
{
	switch (t) {
	
		case NODE_NUMBER:			return "NODE_NUMBER";
		case NODE_VARIABLE:         return "NODE_VARIABLE";
		case NODE_OPERATOR:         return "NODE_OPERATOR";
		case NODE_FUNCTION:         return "NODE_FUNCTION";
		default:					return "UNKNOWN";
	}
}

String print_token(Token* t)
{
	if (is_number(t)) return std::to_string((int)t->num);
	else return t->value;
}

std::string print_ast(AstNode* node) {
    if (node == nullptr) return "";

    std::string result;

    switch (node->type) {
        case NODE_NUMBER:
            result += "NUMBER(" + std::to_string(node->number_value) + ")\n";
            break;

        case NODE_VARIABLE:
            result += "VARIABLE(" + to_string(node->type) + ")\n";
            break;

        case NODE_OPERATOR:
            result += "OPERATOR(" + to_string(node->OperatorNode.op) + ")\n";
            result += "Left: " + print_ast(node->OperatorNode.left);
            result += "Right: " + print_ast(node->OperatorNode.right);
            break;

        default:
            result += "UNKNOWN NODE\n";
            break;
    }

    return result;
}
void print_ast_tree(const AstNode* node, int depth, String prefix, bool first_flag) {
    if (!node) return;

    std::string indent(depth * 4, ' ');  // Create indentation based on depth

    // Print the current node
    switch (node->type) {
        case NODE_NUMBER:
            std::cout << indent << prefix << node->number_value << " NUMBER\n";
            break;
        case NODE_VARIABLE:
            std::cout << indent << prefix << node->symbol << " VARIABLE\n";
            break;
        case NODE_OPERATOR:
            std::cout << indent << prefix << node->symbol << " OPERATOR\n";
            if (node->OperatorNode.left || node->OperatorNode.right) {
				if (first_flag) {
					print_ast_tree(node->OperatorNode.left, depth, "|-- ", false);
					print_ast_tree(node->OperatorNode.right, depth , "`-- ", false);
				} else {
					print_ast_tree(node->OperatorNode.left, depth + 1, "|-- ", false);
					print_ast_tree(node->OperatorNode.right, depth + 1, "`-- ", false);
				}
            }
            break;
        case NODE_FUNCTION:
            std::cout << indent << prefix << node->symbol << " FUNCTION\n";
            print_ast_tree(node->FunctionNode.argument, depth + 1, "`-- ", false);
            break;
        default:
            std::cout << indent << prefix << " UNKNOWN NODE\n";
            break;
    }
}

void print_queue(std::queue<Token*> q)
{
  while (!q.empty()) {
    std::cout << print_token(q.front()) << " ";
    q.pop();
  }
  std::cout << std::endl;
}
//
// DEBUG END
//


