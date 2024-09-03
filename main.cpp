#include <iostream>
#include <cmath>
#include <variant>
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <stack>
#include <queue>
#include <cassert>
#include "glfw3.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944
#define E 2.7182818284590452353602874713527

enum TokenType
{
	NUMBER,
	NUMBER_PI,
	NUMBER_E,

	VARIABLE_X,
	VARIABLE_Y,
	VARIABLE_Z,
	
	PARENTHESIS_OPEN,
	PARENTHESIS_CLOSE,
	OPERATOR_ADD,
	OPERATOR_SUBTRACT,
	OPERATOR_MULTIPLY,
	OPERATOR_DIVIDE,
	OPERATOR_POWER,

	// FUNCTIONS
	FUNCTION_SIN,
	FUNCTION_COS,
	FUNCTION_TAN,
	FUNCTION_CSC,
	FUNCTION_SEC,
	FUNCTION_COT,
	FUNCTION_LOG
};

using String = std::string;

struct Token
{
	String value; // Value associated with the token
	TokenType type; // Type of token

	Token(String v, TokenType t) : value(v), type(t) {}
	Token(double v, TokenType t) : value(std::to_string(v)), type(t) {}
};

std::ostream& operator<<(std::ostream& os, const Token& token) {
	os << token.value;
	return os;
}

std::unordered_map<char, TokenType> char_token_map = {
	{'(', TokenType::PARENTHESIS_OPEN},
	{')', TokenType::PARENTHESIS_CLOSE},

	{'e', TokenType::NUMBER_E},

	// Operators
	{'+', TokenType::OPERATOR_ADD},
	{'-', TokenType::OPERATOR_SUBTRACT},
	{'*', TokenType::OPERATOR_MULTIPLY},
	{'/', TokenType::OPERATOR_DIVIDE},
	{'^', TokenType::OPERATOR_POWER},

	// Valid variable names
	{'x', TokenType::VARIABLE_X},
	{'y', TokenType::VARIABLE_Y},
	{'z', TokenType::VARIABLE_Z}
};

// Quick fix for the implicit mult in lexer. -sep 2 2024
std::unordered_map<char, TokenType> t_map = {
	{'x', TokenType::VARIABLE_X},
	{'y', TokenType::VARIABLE_Y},
	{'z', TokenType::VARIABLE_Z}
};

std::unordered_map<String, TokenType> string_token_map = {
	{"sin", TokenType::FUNCTION_SIN},
	{"cos", TokenType::FUNCTION_COS},
	{"tan", TokenType::FUNCTION_TAN},
	{"csc", TokenType::FUNCTION_CSC},
	{"sec", TokenType::FUNCTION_SEC},
	{"cot", TokenType::FUNCTION_COT},
	{"log", TokenType::FUNCTION_LOG},

	{"pi", TokenType::NUMBER_PI}
};


std::vector<Token> lexer_main(const String &input)
{
	std::vector<Token> tokens;
	size_t i = 0;
	String buffer;
	bool mult_flag = false;

	while (i < input.length()) {
		char c = input[i];
		auto c_it = char_token_map.find(c);

		if (isdigit(c) || c == '.') {
			// Number
			if (isdigit(input[i-1]) || input[i-1] == '.') {
				// Handle multiple digits
				tokens.back().value.append(String(1, c));
			} else {
				// Pushes number token
				tokens.push_back({String(1,c), TokenType::NUMBER});
			}
		} else if (c_it != char_token_map.end()) {
			// variable (x,y,z)
			// Handle multiply without specifying (ie. 3x == 3*x)
			if (isdigit(input[i-1]) && t_map.find(c) != t_map.end()) {
				tokens.push_back({String(1,'*'), TokenType::OPERATOR_MULTIPLY});
			}
			tokens.push_back({String(1,c), c_it->second});

		} else if (isalpha(c)) {
			// function
			buffer.push_back(c);
			auto s_it = string_token_map.find(buffer);
			if (isdigit(input[i-1]) && buffer.size() == 1) tokens.push_back({String(1,'*'), TokenType::OPERATOR_MULTIPLY});
			if (s_it != string_token_map.end()) {
				tokens.push_back({buffer, s_it->second});
				buffer.clear();
			}

		} else {
			buffer.clear();
		}

		i++;
	}
	return tokens;
}

using BinaryFunction = double(*)(double, double);
using UnaryFunction = double(*)(double);

// MATH FUNCTIONS
// BASIC
double power(double a, double b) { return std::pow(a, b); }
double multiply(double a, double b) { return a * b; }
double divide(double a, double b) { if (b == 0) throw std::runtime_error("Division by zero"); return a / b; }
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
double math_log(double a) { if (a <= 0) throw std::runtime_error("Log of zero or negative number"); return std::log(a); }

struct Operator
{
	TokenType type;
	int precedence;
	bool associativity; // true is left
	BinaryFunction f;
};

std::unordered_map<TokenType, Operator> operator_map = {
	{TokenType::OPERATOR_POWER, {TokenType::OPERATOR_POWER, 4, false, power}},

	{TokenType::OPERATOR_MULTIPLY, {TokenType::OPERATOR_MULTIPLY, 3, true, multiply}},
	{TokenType::OPERATOR_DIVIDE, {TokenType::OPERATOR_DIVIDE, 3, true, divide}},

	{TokenType::OPERATOR_ADD, {TokenType::OPERATOR_ADD, 2, true, add}},
	{TokenType::OPERATOR_SUBTRACT, {TokenType::OPERATOR_SUBTRACT, 2, true, subtract}}

};

std::unordered_map<TokenType, UnaryFunction> unary_function_map = {
	{TokenType::FUNCTION_SIN, math_sin},
	{TokenType::FUNCTION_COS, math_cos},
	{TokenType::FUNCTION_TAN, math_tan},
	{TokenType::FUNCTION_CSC, math_csc},
	{TokenType::FUNCTION_SEC, math_sec},
	{TokenType::FUNCTION_COT, math_cot},
	{TokenType::FUNCTION_LOG, math_log}
};



bool is_function(Token t) { return (t.type >= FUNCTION_SIN && t.type <= FUNCTION_LOG); }
bool is_variable(Token t) { return (t.type >= VARIABLE_X && t.type <= VARIABLE_Z); }
bool is_operator(Token t) { return (t.type >= OPERATOR_ADD && t.type <= OPERATOR_POWER); }
bool is_number(Token t) { return (t.type >= NUMBER && t.type <= NUMBER_E); }

std::queue<Token> parser_main(const std::vector<Token> &tokens)
{

	std::queue<Token> output;
	std::stack<Token> operators;


	// Shunting yard algorithm, see https://en.wikipedia.org/wiki/Shunting_yard_algorithm
	for (const Token &token : tokens) {

		// Numbers go straight to output
		if (is_number(token)) {
			if (token.type == NUMBER_PI) {
				output.push({PI, token.type});
			} else if (token.type == NUMBER_E) {
				output.push({E, token.type});
			} else {
				output.push(token);
			}
		} else if (is_variable(token)) {
			output.push(token);
		} else if (is_function(token)) {
			operators.push(token);
		} else if (is_operator(token)) {
			auto o1_it = operator_map.find(token.type);
			Operator o_1 = o1_it->second;

			while (!operators.empty() && is_operator(operators.top()) && operators.top().type != TokenType::PARENTHESIS_OPEN) {
				auto o2_it = operator_map.find(operators.top().type);
				Operator o_2 = o2_it->second;
				if (o_2.precedence > o_1.precedence || (o_1.precedence == o_2.precedence && o_1.associativity)) {
					output.push(operators.top());
					operators.pop();
				} else break;
			}
			operators.push(token);
		} else if (token.type == TokenType::PARENTHESIS_OPEN) {
			operators.push(token);
		} else if (token.type == TokenType::PARENTHESIS_CLOSE) {
			while (!operators.empty() && operators.top().type != TokenType::PARENTHESIS_OPEN) {
				assert(!operators.empty()); // missing matching parenthesis
				output.push(operators.top());
				operators.pop();
			}
			assert(operators.top().type == TokenType::PARENTHESIS_OPEN); // missing matching parenthesis
			operators.pop();
			if (!operators.empty() && is_function(operators.top())) {
				output.push(operators.top());
				operators.pop();
			}
		}
	}

	while (!operators.empty()) {
		assert(operators.top().type != TokenType::PARENTHESIS_OPEN); // missing matching parenthesis
		output.push(operators.top());
		operators.pop();
	}

	return output;
}

std::unordered_map<char, double> variable_map = {
	{'x', 1.0},
	{'y', 1.0},
	{'z', 1.0}
};


// No variable
double operate(const std::vector<Token> &t, const Token& f)
{
	if (is_operator(f)) {
		assert(t.size() == 2);
		auto it = operator_map.find(f.type);
		Operator o = it->second;
		double result;
		try {
            double left_operand = std::stod(t[0].value);
            double right_operand = std::stod(t[1].value);
			// @Note : Associativity is dealt with in the algorithm. -sep 2 2024
			result = o.f(left_operand, right_operand);

        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for std::stod: " << e.what() << " " << t[0].value << " " << t[1].value << std::endl;
            throw;
        }



		return result;

	} else if (is_function(f)) {
		assert(t.size() == 1);

		auto it = unary_function_map.find(f.type);
		if (it != unary_function_map.end()) {
			auto func = it->second;
			double input = std::stod(t[0].value);
            double result = func(input);
            return result;
			
		}
	} else throw std::runtime_error("not supported yet");
	return 0.0;
}

double solve_main(std::queue<Token> tokens)
{
	std::stack<Token> temp;

	while (!tokens.empty()) {
		const auto& token = tokens.front();

		if (is_number(token)) {
			temp.push(token);
			tokens.pop();
			continue;
		} else if (is_variable(token)) {
			temp.push({variable_map.at(token.value[0]), token.type});
			tokens.pop();
			continue;
		}

		std::vector<Token> operands;
		if (is_operator(token)) {
			assert(temp.size() >= 2);
			Token right_operand = temp.top();
			temp.pop();
			Token left_operand = temp.top();
			temp.pop();
			operands.push_back(left_operand);
			operands.push_back(right_operand);
		} else if (is_function(tokens.front())) {
			assert(temp.size() >= 1);
			operands.push_back(temp.top());
			temp.pop();
		} else throw std::runtime_error("invald token encountered.");

		double result = operate(operands, token);
		tokens.pop();
		temp.push({std::to_string(result), TokenType::NUMBER});
	}
	// Sanity check
	assert(temp.size() == 1);
	return std::stod(temp.top().value);
}

std::vector<std::pair<double, double>> evaluate_function(const std::queue<Token>& tokens, double start, double end, double step)
{
	std::vector<std::pair<double, double>> results;
	// Iterate through the range of values
	for (double x = start; x <= end; x+=step) {
		// Update map
		variable_map['x'] = x;

		// Evaluate function
		double y = solve_main(tokens);

		results.emplace_back(x, y);

	}

	return results;
}

void plot_function(const std::queue<Token>& tokens, double start, double end, double step)
{
	
}

//@Debug
void print_queue(const std::queue<Token> &q)
{
	std::queue<Token> temp = q;
	while (!temp.empty()) {
		const Token &token = temp.front();
		std::cout << token << " ";
		temp.pop();
	}
	std::cout << "\n";
}

void framebugger_size_callback(GLFWwindow*, int width, int height)
{
	glViewport(0,0,width,height);
}

int main()
{
	//String test_string = "(sin(x)+x^2*cos(x)) / (e^x-log(x))";
	String test_string = "(sin(x)+x^2*cos(x)) / (e^x-log(x))";
	std::vector<Token> tokens = lexer_main(test_string);
	// @Debug
	std::cout << "Input-:\n" << test_string << "\n\n";
	std::cout << "Lexer-:\n";
	for (Token t : tokens) std::cout << t << " " << t.type << "\n";
	std::cout << "\nParser:-\n";
	for (Token t : tokens) std::cout << t;
	std::cout << "\n";
	std::queue<Token> parsed = parser_main(tokens);
	std::cout << "\nAlgorithm:-\n";
	print_queue(parsed);
	std::cout << "\n\n";
	std::cout << solve_main(parsed) << std::endl;
	std::cout << "Plot:-\n";
	evaluate_function(parsed, 0.01, 2, 0.01);
	return 0;
}
