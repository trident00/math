#include "math.h"
#include <initializer_list>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include <thread>
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
#include "raylib.h"
#include <thread>
#include <limits>
#include "parser.h"
#include "ast.h"

#define E		2.7182818284590452353602874713527
#define INF		std::numeric_limits<float>::max()
#define FINISH() while (!WindowShouldClose()) { std::cout<<"sdfsdsdf"<<std::endl; }
#define For(container) for (auto& item : container)


std::unordered_map<char, double> variable_map = {
	{'x', 1.0f},
	{'y', 1.0f},
	{'z', 1.0f}
};

typedef std::pair<bool, int> PlotState;

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
	std::vector<std::vector<Vector2>> grid_coordinates;
	double xmi, xmx, dx; // Config
    double x_min, y_min, x_max, y_max;
    double x_scale, y_scale, x_offset, y_offset;
	Vector2 max_coords = {graph_width, graph_height};
	float axis_limit = 99999999.0f;
};

std::vector<Vector2> evaluate_function(AstNode* tree, double start, double end, double step)
{
	std::vector<Vector2> results;
	// Iterate through the range of values
	for (double x = start; x <= end+step; x+=step) {
		// Update map
		variable_map['x'] = x;
		// Evaluate function
		double y = ast_solve(tree, x);
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
std::vector<Vector2> cal_axes(PlotData& data, float x, float y)
{
	std::vector<Vector2> a = {{(float)data.xmi, y}, {(float)data.xmx, y}, {x, -data.axis_limit}, {x, data.axis_limit}};
	return standard_coords(data, a);
}
void cal_grid(PlotData& data, int grid_step)
{
	// X
	std::vector<Vector2> x_points; std::vector<Vector2> y_points;
	
	// Correct height of grid
	For (x_points) {
		printf("(%f, %f)\n", item.x, item.y);
		//printf("(%f, %f)->(%f, %f) : (%f, %f)->(%f, %f)\n", item[0].x, item[0].y, item[1].x, item[1].y, item[2].x, item[2].y, item[3].x, item[3].y);
	}
	For (y_points) {
		printf("(%f, %f)\n", item.x, item.y);
	}
}
//
//	DRAW FUNCTIONS
//
void draw_axes(PlotData& data, const Color& color)
{
	auto vec = cal_axes(data, 0, 0);
	// Correct height of axis
	DrawLineV(vec[0], vec[1], color);
	DrawLineV(vec[2], vec[3], color);
}
void draw_grid(PlotData& data, PlotState& state, const Color& color, int x)
{
		auto vec = data.grid_coordinates[x];
		DrawLineV(vec[0], vec[1], color);
		DrawLineV(vec[2], vec[3], color);
}
void draw_function(PlotData& data, PlotState& state, const Color& color)
{
	if (state.second > 1 && !data.screen_coordinates.empty()) {
		//printf("(%.10f, %.10f)\n",data.screen_coordinates.back().x ,data.screen_coordinates.back()		std::numeric_limits<float>::max().y);
		int num_points = std::min(state.second, static_cast<int>(data.screen_coordinates.size()));
		DrawLineStrip(data.screen_coordinates.data(), num_points, color);
	}
}

void update_drawing(PlotData& data, PlotState& state, int points, int fps, float total_time, bool s)
{
	float points_per_frame = static_cast<float>(points) / (fps * total_time * 30);
	if (state.second < points) {
		state.second += points_per_frame;
		if (state.second > points) {
			state.second = points;
			state.first = true;
		}
	}
}

// One function
void prepare_data(PlotData& data)
{
	cal_minmax(data);
	cal_scale(data);
	cal_coords(data);
	cal_axes(data,0,0);
	cal_grid(data, 1);
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

	std::vector<PlotState> graph_state(plots.size(), {false, 0});
	std::vector<PlotState> fn_state(plots.size(), {false, 0});
	while (!WindowShouldClose()) {
		ClearBackground(Color());
		BeginDrawing();
		for (int x = 0; x < (int)plots[first]->grid_coordinates.size(); x++) {
			draw_grid(*plots[first], graph_state[first], {118,118,118,255}, x);
			//update_drawing(*plots[first], graph_state[first], plots[first]->grid_coordinates.size(), fps, true);
		}
		bool all_functions_drawn = true;
		draw_axes(*plots[first], {214,214,214,255});
		for (size_t i = 0; i < plots.size(); i++) {
			draw_function(*plots[i], fn_state[i], colors[i]);
			update_drawing(*plots[i], fn_state[i], plots[i]->screen_coordinates.size(), fps, 2, false);
			if (!fn_state[first].first) all_functions_drawn = false;
		}

		if (all_functions_drawn) {
			draw_axes(*plots[first], {214,214,214,255});
		}
		EndDrawing();
	}
}

PlotData function(String f, double xmi, double xmx, double dx)
{
	PlotData data;
	data.function = f; data.xmi = xmi; data.xmx = xmx; data.dx = dx;
	data.coordinates = evaluate_function(ast_tree(f), xmi, xmx, dx);
	return data;
}

int main()
{
	InitWindow(1000, 800, "Math.exe");
	PlotData f_2 = function("0.1xx", -4, 4, 0.001);
	PlotData f_1 = function("cos(sin(tanx))", -4, 4, 0.001);
	std::vector<PlotData*> fs = {&f_1, &f_2};
	prepare_data(f_1);
	prepare_data(f_2, f_1);
	plot(fs);
	return 0;
}
