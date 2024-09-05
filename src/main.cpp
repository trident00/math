#include "math.h"
#include <iostream>
#include <chrono>
#include <array>
#include <utility>
#include <cmath>
#include <variant>
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <stack>
#include <queue>
#include <cassert>
#include <GLFW/glfw3.h>
#include <thread>

#define PI 3.141592653589793238462643383279502884197169399375105820974944
#define E 2.7182818284590452353602874713527

#define FINISH() while (!glfwWindowShouldClose(window)) { glfwPollEvents(); }


using String = std::string;

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
double math_log(double a) { if (a <= 0); return std::log(a); }

// COLOR
struct Color
{
	int r;
	int g;
	int b;
};

void color(const Color& color) { glColor3f(color.r / 255.0f, color.g / 255.0f, color.b / 255.0f); }

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


GLFWwindow* init_window()
{
	// Initialize the library
	if (!glfwInit()) throw std::runtime_error("library failed");

	// Set the window to be non-resizable
	glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

	// Create windowed mode window and its OpenGL context
	String title = "Math.exe";
	GLFWwindow* window = glfwCreateWindow(1000, 800, title.c_str(), NULL, NULL);
	if (!window) {
		glfwTerminate();
	}

	// Make window's context current
	glfwMakeContextCurrent(window);
	// Register the framebuffer size callback
	glfwSetFramebufferSizeCallback(window, [](GLFWwindow*, int, int) {
		glViewport(0, 0, 800, 800);
	});

	return window;
}

int end_window()
{
	glfwTerminate();
	return 0;
}

GLFWwindow* window = init_window();

std::pair<double, double> to_ndc(std::pair<double, double> point,
				double x_scale, double y_scale,
				double x_offset, double y_offset)
{
	double x_ndc = point.first * x_scale + x_offset;
	double y_ndc = point.second * y_scale + y_offset;

	return {x_ndc, y_ndc};
}


struct PlotState {
	int current_point = 0;
	bool drawing_complete = false;
};

struct PlotData {
	String function;
	std::vector<std::pair<double, double>> coordinates;
	bool smooth;
	bool incremental;

	// Useful values
	double x_min, y_min, x_max, y_max, x_range, y_range, x_scale, y_scale, x_offset, y_offset;
	
	// Axes coordinates
	std::vector<std::pair<double, double>> axes_coordinates;

	// To plot points to screen
	std::vector<std::pair<double, double>> ndc_coordinates;
	
	PlotData(const String& func, const auto& coords, bool sm, bool inc)
	: function(func), coordinates(coords), smooth(sm), incremental(inc) {}

	//
	// PREPARE DATA
	//
	void c_minmax()
	{
		// Min/Max values
		x_min = std::numeric_limits<double>::max();
		x_max = std::numeric_limits<double>::lowest();
		y_min = std::numeric_limits<double>::max();
		y_max = std::numeric_limits<double>::lowest();

		for (const auto& point : coordinates) {
			x_min = std::min(x_min, point.first);
			x_max = std::max(x_max, point.first);
			y_min = std::min(y_min, point.second);
			y_max = std::max(y_max, point.second);
		}
	}
	void c_range()
	{
		// Mapping to NDC for use of GLFW3
		x_range = x_max - x_min;
		y_range = y_max - y_min;
	}
	void c_scale()
	{

		// NOTE: MAKE THE SCALE LIMITED AT SOME POINT, AS TO VISUALIZE ASYMPTOTIC/EXPONENTIONAL FUNCTIONS
		// Calculate adjusted ranges
		double scale = std::max(x_range, y_range);
		x_scale = 2.0 / scale;
		y_scale = 2.0 / scale;

		x_offset = -(x_min + x_range / 2.0) * x_scale;
		y_offset = -(y_min + y_range / 2.0) * y_scale;
	}
	void c_first() { c_minmax(); c_range(); c_scale(); }
	void c_next() { c_minmax(); c_range(); }

	//
	// COORDINATES
	//
	void c_axes()
	{
		axes_coordinates = {
			{to_ndc({x_min, 0}, x_scale, y_scale, x_offset, y_offset).first,	
			to_ndc({x_min, 0}, x_scale, y_scale, x_offset, y_offset).second},
			{to_ndc({x_max, 0}, x_scale, y_scale, x_offset, y_offset).first,
			to_ndc({x_max, 0}, x_scale, y_scale, x_offset, y_offset).second},		
			{to_ndc({0, y_min}, x_scale, y_scale, x_offset, y_offset).first,
			to_ndc({0, y_min}, x_scale, y_scale, x_offset, y_offset).second},
			{to_ndc({0, y_max}, x_scale, y_scale, x_offset, y_offset).first,
			to_ndc({0, y_max}, x_scale, y_scale, x_offset, y_offset).second}
		};
	}
	void c_ndc()
	{
		for (const auto& point : coordinates) {
			ndc_coordinates.emplace_back(to_ndc(point, x_scale, y_scale, x_offset, y_offset));
		}
	}

	//
	// DRAW
	//
	void draw_axes(const Color& rgb) const
	{
		color(rgb);
		glBegin(GL_LINES);
		glVertex2f(axes_coordinates[0].first, axes_coordinates[0].second);
		glVertex2f(axes_coordinates[1].first, axes_coordinates[1].second);
		glVertex2f(axes_coordinates[2].first, axes_coordinates[2].second);
		glVertex2f(axes_coordinates[3].first, axes_coordinates[3].second);
		glEnd();
	}
	void draw_function(const Color& rgb, PlotState& state) const
	{
		const int seconds = 2;
		int points_per_frame = ndc_coordinates.size() / (seconds * 120);
		
		color(rgb);
		glBegin(smooth ? GL_LINE_STRIP : GL_POINTS);
		for (int i = 0; i < state.current_point && i < ndc_coordinates.size(); ++i) {
			glVertex2f(ndc_coordinates[i].first, ndc_coordinates[i].second);
		}

		glEnd();
		// Check if drawing is complete
		if (state.current_point < ndc_coordinates.size()) {
			state.current_point += points_per_frame;
		} else {
			state.drawing_complete = true;
		}
	}
};

void prepare_data(PlotData& plot_data, bool main)
{
	if (main) {
		plot_data.c_first();
		plot_data.c_axes();
		plot_data.c_ndc();
	}
	else {
		plot_data.c_next();
	}
}

auto plot_function(const PlotData& plot_data, bool other) 
{
	PlotState plot_state;
	glViewport(0,0,800,800);
	while (!glfwWindowShouldClose(window)) {
		// Axes
		plot_data.draw_axes({214,214,214});
		// Points
		plot_data.draw_function({255,0,0}, plot_state);
		// Swap front and back buffers
		glfwSwapBuffers(window);
		// Poll for and process events
		glfwPollEvents();

		if (plot_state.drawing_complete) {
			// Keep window open until closed
			FINISH();
		}
	}
}

void draw_points(PlotData& plot_data, const Color& rgb)
{
	color(rgb);
	
}

void plot_functions(PlotData& plot_data1, const PlotData& plot_data2)
{
	plot_function(plot_data1, true);
}

void menu()
{
	glViewport(800, 0, 200, 800);
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

void plot(String function)
{
	const auto& parsed = parser_main(lexer_main(function));
	const auto& points = evaluate_function(parsed, -4, 4, 0.01);
	PlotData plot_data = {function, points, true, true};
	prepare_data(plot_data, true);
	plot_function(plot_data, false);
}

int main()
{
	//String test_string = "(sin(x)+x^2*cos(x)) / (e^x-log(x))";
	//String test_string = "sin(x)";
	//std::vector<Token> tokens = lexer_main(test_string);
	//// @Debug
	//std::cout << "Input-:\n" << test_string << "\n\n";
	//std::cout << "Lexer-:\n";
	//for (Token t : tokens) std::cout << t << " " << t.type << "\n";
	//std::cout << "\nParser:-\n";
	//for (Token t : tokens) std::cout << t;
	//std::cout << "\n";
	//std::cout << "\nAlgorithm:-\n";
	//print_queue(parsed);
	//std::cout << "\n\n";
	//std::queue<Token> parsed = parser_main(tokens);
	//const auto& points = evaluate_function(parsed, -2, 2, 0.01);
	//PlotData plot_data(test_string, points, true, true);
	//std::cout << "# of coords:-\n" << points.size();
	plot("sin(x^2)*cos(x^2)/sin(x)");
	end_window();
	return 0;
}
