#include "math.h"
#include <initializer_list>
#include <iostream>
#include <stdio.h>
#include <chrono>
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
//#include <GLFW/glfw3.h>
#include "raylib.h"
#include <thread>

// PI defined in raylib
//#define PI 3.141592653589793238462643383279502884197169399375105820974944
#define E 2.7182818284590452353602874713527

#define FINISH() while (!WindowShouldClose()) { std::cout<<"sdfsdsdf"<<std::endl; }
//#define For(container) for (auto& item : container)

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

std::vector<Vector2> evaluate_function(const std::queue<Token>& tokens, double start, double end, double step)
{
	std::vector<Vector2> results;
	// Iterate through the range of values
	for (double x = start; x <= end+step; x+=step) {
		// Update map
		variable_map['x'] = x;
		// Evaluate function
		double y = solve_main(tokens);
		//printf("(%.10f, %.10f)\n", x, y);
		Vector2 vec = {(float)x,(float)y};
		results.emplace_back(vec);
	}
	return results;
}

const int screen_width	= 1000;
const int screen_height	= 800;
const int graph_width	= 800;
const int graph_height	= 800;

Vector2 screen_coords(Vector2 vec, double x_scale, double y_scale, double x_offset, double y_offset)
{
	double x = (vec.x * x_scale) + x_offset;
	double y = (vec.y * y_scale) + y_offset;
	// converts NDC to regular coords (laziness)
	double standard_x = graph_width * (0.5*x + 0.5);
	double standard_y = graph_width * (-0.5*y + 0.5);
	return {(float)standard_x, (float)standard_y};
}

struct PlotState {
};

struct PlotData {
	String function;
	std::vector<Vector2> coordinates;
	bool smooth;
	bool incremental;

	int current_point = 0;
	bool drawing_complete = false;
	// Useful values
	double x_min, y_min, x_max, y_max, x_range, y_range, x_scale, y_scale, x_offset, y_offset;

	// To plot points to screen
	std::vector<Vector2> screen_coordinates, axes_coordinates;
	
	PlotData(const String& func, const std::vector<Vector2>& coords, bool sm, bool inc)
	: function(func), coordinates(coords), smooth(sm), incremental(inc) { coordinates = coords; }

	//
	// PREPARE DATA
	//
	void c_minmax()
	{
		// Min/Max values
		x_min = std::numeric_limits<float>::max();
		x_max = std::numeric_limits<float>::lowest();
		y_min = std::numeric_limits<float>::max();
		y_max = std::numeric_limits<float>::lowest();

		for (const auto& vec : coordinates) {
			x_min = std::min((float)x_min, vec.x);
			x_max = std::max((float)x_max, vec.x);
			y_min = std::min((float)y_min, vec.y);
			y_max = std::max((float)y_max, vec.y);
		}
		//printf("Init: %f : %f : %f : %f : %f : %f : %f : %f\n\n",x_min,x_max,y_min,y_max, x_scale, y_scale, x_offset, y_offset);
	}
	void c_range()
	{
		// Mapping to NDC for use of GLFW3
		x_range = x_max - x_min;
		y_range = y_max - y_min;
	}
	void c_scale()
	{
		// NOTE: MAKE THE SCALE LIMITED AT SOME POINT, AS TO VISUALIZE ASYMPTOTIC/EXPONENTIAL FUNCTIONS

		// Normalization of coords (shift and scale)
		double scale = std::max(x_range, y_range);
		x_scale = 2.0 / scale;
		y_scale = 2.0 / scale;
		
		x_offset = -(x_min + x_range / 2.0) * x_scale;
		y_offset = -(y_min + y_range / 2.0) * y_scale;
	}
	void c_first() { c_minmax(); c_range(); c_scale();}

	void c_next(PlotData* main) {
		c_minmax(); c_range();
		x_scale = main->x_scale;
		y_scale = main->y_scale;
		x_offset = main->x_offset;
		y_offset = main->y_offset;
	}

	//
	// COORDINATES
	//
	void c_axes()
	{
		axes_coordinates = {
			{screen_coords({(float)x_min, 0}, x_scale, y_scale, x_offset, y_offset)},	
			{screen_coords({(float)x_max, 0}, x_scale, y_scale, x_offset, y_offset)},
			{screen_coords({0, -10000}, x_scale, y_scale, x_offset, y_offset)},
			{screen_coords({0,  10000}, x_scale, y_scale, x_offset, y_offset)}
		};
	}

	void c_ndc()
	{
		//printf("Init: %f : %f : %f : %f : %f : %f : %f : %f\n\n",x_min,x_max,y_min,y_max, x_scale, y_scale, x_offset, y_offset);
		for (const auto& point : coordinates) {
			//printf("(%.10f, %.10f) : ", point.x, point.y);
			screen_coordinates.emplace_back(screen_coords(point, x_scale, y_scale, x_offset, y_offset));
			//printf("(%.10f, %.10f)\n",screen_coordinates.back().x ,screen_coordinates.back().y);
		}
		//printf("Post: %f : %f : %f : %f : %f : %f : %f : %f\n\n",x_min,x_max,y_min,y_max, x_scale, y_scale, x_offset, y_offset);
	}

	//
	// DRAW
	//
	void draw_axes(const Color& color) const
	{
		DrawLineV(axes_coordinates[0], axes_coordinates[1], color);
		DrawLineV(axes_coordinates[2], axes_coordinates[3], color);
	}
	void draw_function(const Color& color)
	{
		const int seconds = 2;
		int points_per_frame = (screen_coordinates.size() / (seconds * 60));

		if (current_point > 1 && !screen_coordinates.empty()) {
			int num_points = std::min(current_point, static_cast<int>(screen_coordinates.size()));

			// Filter coordinates if needed.	
			std::vector<Vector2> filtered_coordinates;
			for (int i = 0; i < num_points; i++) {
				// keep x coord from going beyond the boundary of the graph (x_max, inf)
				if (!(screen_coordinates[i].x > graph_width)) {
					filtered_coordinates.push_back(screen_coordinates[i]);
				}
			}
	
			DrawLineStrip(filtered_coordinates.data(), filtered_coordinates.size(), color);
		}

		// Check if drawing is complete
		if (current_point < screen_coordinates.size()) {
			current_point += points_per_frame;
		} else {
			drawing_complete = true;
		}
	}
};

// One function
void prepare_data(PlotData* plot_data)
{
	plot_data->c_first();
	plot_data->c_axes();
	plot_data->c_ndc();
}

// Two functions, secondary function takes the scale of the main
void prepare_data(PlotData* plot_data, PlotData* plot_data_main)
{
	plot_data->c_next(plot_data_main);
	plot_data->c_ndc();
}

struct Function
{
	String str;
	double x_min, x_max, delta_x;
	PlotData* plot_data;

	Function(const String& s, const double& xmi, const double& xmx, const double& dx)
	: str(s), x_min(xmi), x_max(xmx), delta_x(dx), plot_data(nullptr)
	{
		auto points = evaluate_function(parser_main(lexer_main(str)), x_min, x_max, delta_x);
		plot_data = new PlotData(str, points, true, true);
	}

	void clean()
	{
		delete plot_data;
		plot_data = nullptr;
	}
};

void plot(std::vector<Function>& functions, int first)
{
	const std::vector<Color> colors = {Color(255,0,0,255), Color(0,255,0,255), Color(0,0,255,255), Color(255,0,255,255)};

	// Calculate data
	for (size_t i = 0; i < functions.size(); i++) {
		// Create PlotData object
		auto* plot = functions[i].plot_data;
		// Prepare function data
		if (i == first) {
			prepare_data(plot);
		}
		else { 
			prepare_data(functions[i].plot_data, functions[first].plot_data); 
		}
	}

	// Draw
	SetTargetFPS(60);
	while (!WindowShouldClose()) {
		ClearBackground(Color());
		BeginDrawing();
		// START
		bool all_functions_drawn = true;
		for (size_t i = 0; i < functions.size(); i++) {
			//if (i == first) plot->draw_axes({214,214,214,255});
			functions[i].plot_data->draw_function(colors[i]);

			if (!functions[i].plot_data->drawing_complete) all_functions_drawn = false;
		}
		// After Functions drawn
		if (all_functions_drawn) {
			functions[first].plot_data->draw_axes({214,214,214,255});
		}
		// END
		EndDrawing();
	}
}

int main()
{
	Function f_1 = {"e^(sinx)", -4, 4, 0.001};
	Function f_2 = {"1-x^2/2-x^4/8+3*x^6/80", -4, 4, 0.001};
	Function f_3 = {"4-0.04x^2", -4, 4, 0.001};
	Function f_4 = {"cos(sin(tanx))", -4, 4, 0.0001};

	InitWindow(1000, 800, "Math.exe");
	std::vector<Function> fs = {f_4, f_2, f_1, f_3};
	std::cout << std::endl;
	plot(fs, 0);
	f_1.clean(); f_2.clean(); f_3.clean(); f_4.clean();
	//plot("cos(sin(tanx))", "1-x^2/2-x^4/8+(3*x^6)/80");
	//plot("cos(sin(tanx))", "1-x^2/2-x^4/8");
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

	return 0;
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
