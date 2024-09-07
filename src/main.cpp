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
#define For(container) for (auto& item : container)

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

const int screen_width	= 1000;
const int screen_height	= 800;
const int graph_width	= 800;
const int graph_height	= 800;
// Types and Data Structures
struct PlotData {
    std::string function;
    std::vector<Vector2> coordinates;
    std::vector<Vector2> axes_coordinates;
    std::vector<Vector2> screen_coordinates;
	double xmi, xmx, dx; // Config
    double x_min, y_min, x_max, y_max;
    double x_scale, y_scale, x_offset, y_offset;
    int current_point = 0;
    bool drawing_complete = false;
};

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

//
//	CALCULATE FUNCTIONS
//
void cal_minmax(PlotData& data)
{
	data.x_min = std::numeric_limits<float>::max();
	data.x_max = std::numeric_limits<float>::lowest();
	data.y_min = std::numeric_limits<float>::max();
	data.y_max = std::numeric_limits<float>::lowest();

	for (const auto& vec : data.coordinates) {
		data.x_min = std::min((float)data.x_min, vec.x);
		data.x_max = std::max((float)data.x_max, vec.x);
		data.y_min = std::min((float)data.y_min, vec.y);
		data.y_max = std::max((float)data.y_max, vec.y);
	}
}
void cal_scale(PlotData& data)
{
	double x_range = data.x_max - data.x_min;
	double y_range = data.y_max - data.y_min;
	// Normalization of coords
	double scale = std::max(x_range, y_range);
	data.x_scale = 2.0 / scale;
	data.y_scale = 2.0 / scale;

	data.x_offset = -(data.x_min + x_range / 2.0) * data.x_scale;
	data.y_offset = -(data.y_min + y_range / 2.0) * data.y_scale;
}

Vector2 to_coords(PlotData& data, Vector2& point)
{
	double ndc_x = (point.x * data.x_scale) + data.x_offset;
	double ndc_y = (point.y * data.y_scale) + data.y_offset;
	// converts NDC to regular coords (laziness)
	double standard_x = graph_width * (0.5*ndc_x + 0.5);
	double standard_y = graph_height * (-0.5*ndc_y + 0.5);
	return {float(standard_x), float(standard_y)};
}

std::vector<Vector2> standard_coords(PlotData& data, std::vector<Vector2>& points)
{
	std::vector<Vector2> result;
	for (auto& point : points) {
		result.emplace_back(to_coords(data, point));
	}
	return result;
}

void cal_coords(PlotData& data)
{
	auto vec = standard_coords(data, data.coordinates);
	For (vec) {
		// TEST CASES
		// keep x coord from going beyond the boundary of the graph (x_max, inf)
		if (item.x <= graph_width) {
			data.screen_coordinates.emplace_back(item);
		}
	}
}
void cal_axes(PlotData& data)
{
	std::vector<Vector2> a = {{(float)data.x_min, 0}, {(float)data.x_max, 0}, {0, -10000}, {0, 10000}};
	data.axes_coordinates = standard_coords(data, a);
	std::cout<<data.axes_coordinates.size()<<std::endl;
}
//
//	DRAW FUNCTIONS
//
void draw_axes(const PlotData& data, const Color& color)
{
	DrawLineV(data.axes_coordinates[0], data.axes_coordinates[1], color);
	DrawLineV(data.axes_coordinates[2], data.axes_coordinates[3], color);
}
void draw_function(PlotData& data, const Color& color)
{
	if (data.current_point > 1 && !data.screen_coordinates.empty()) {
		//printf("(%.10f, %.10f)\n",data.screen_coordinates.back().x ,data.screen_coordinates.back().y);
		int num_points = std::min(data.current_point, static_cast<int>(data.screen_coordinates.size()));
		DrawLineStrip(data.screen_coordinates.data(), num_points, color);
	}
}

void update_drawing(PlotData& data, int fps=60)
{
	int points_per_frame = (data.screen_coordinates.size() / (2 * fps));
	if (data.current_point < (int)data.screen_coordinates.size()) {
		data.current_point += points_per_frame;
	} else {
		data.drawing_complete = true;
	}
}

// One function
void prepare_data(PlotData& data)
{
	cal_minmax(data);
	cal_scale(data);
	cal_coords(data);
	cal_axes(data);
}
// Two functions, secondary function takes the scale of the main
void prepare_data(PlotData& data, PlotData& data_main)
{
	cal_minmax(data);
	data.x_scale = data_main.x_scale;
	data.y_scale = data_main.y_scale;
	data.x_offset = data_main.x_offset;
	data.y_offset = data_main.y_offset;
	cal_coords(data);
}

void plot(std::vector<PlotData*>& plots, int first=0, int fps=60)
{

	const std::vector<Color> colors = {Color(255,0,0,255), Color(0,255,0,255), Color(0,0,255,255), Color(255,0,255,255)};

	std::cout << plots[0]->screen_coordinates.size() << std::endl;
	std::cout << plots[1]->screen_coordinates.size() << std::endl;
//	for (const auto& coord : plots[0].screen_coordinates) {
//		std::cout << "Screen Coord: (" << coord.x << ", " << coord.y << ")" << std::endl;
//	}
	while (!WindowShouldClose()) {
		ClearBackground(Color());
		BeginDrawing();
		bool all_functions_drawn = true;
		for (size_t i = 0; i < plots.size(); i++) {
			draw_function(*plots[i], colors[i]);
			update_drawing(*plots[i]);

			if (!plots[i]->drawing_complete) all_functions_drawn = false;
			//plots[i]->draw_axes({214,214,214,255});
		}

		if (all_functions_drawn) {
			draw_axes(*plots[first],  {214,214,214,255});
		}
		EndDrawing();
	}
}

PlotData function(String f, double xmi, double xmx, double dx)
{
	PlotData data;
	data.function = f; data.xmi = xmi; data.xmx = xmx; data.dx = dx;
	data.coordinates = evaluate_function(parser_main(lexer_main(f)), xmi, xmx, dx);
	return data;
}

int main()
{
	InitWindow(1000, 800, "Math.exe");
	PlotData f_1 = function("cos(sin(tanx))", -4, 4, 0.0001);
	PlotData f_2 = function("1-x^2/2-x^4/8+3*x^6/80", -4, 4, 0.001);
	std::vector<PlotData*> fs = {&f_1, &f_2};
	prepare_data(f_1);
	prepare_data(f_2, f_1);
	std::cout << f_1.screen_coordinates.size() << std::endl;
	std::cout << f_2.screen_coordinates.size() << std::endl;
	plot(fs);
	return 0;
}
