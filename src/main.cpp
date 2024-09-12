#include "parser.h"
#include "raylib.h"
#include "ast.h"
#include <iostream>
#include <stdio.h>
#include <cassert>
#include <limits>

#define E		2.7182818284590452353602874713527
#define INF     std::numeric_limits<double>::max()

const int graph_width = 1000;
const int graph_height = 1000;

struct Coordinates2D
{
	std::vector<Vector2> xy_coordinates;
	std::pair<double, double> domain;
	std::pair<double, double> range;
	std::vector<Vector2> screen_coordinates;
	std::vector<Vector2> sa_coords; // Screen Axes Coords
	Coordinates2D() {}
	~Coordinates2D() {}
};

struct Function
{
	String str;
	AstNode* function; // ast_tree(function)
	Coordinates2D* coords;
	
	Function() {}
	~Function() {}
};

Vector2 normalize(double scale, double x, double y, double x_offset, double y_offset)
{
	double ndc_x = (x * scale) + x_offset;
	double ndc_y = (y * scale) + y_offset;
	// converts NDC to regular coords (laziness)
	double screen_x = graph_width * (0.5*ndc_x + 0.5);
	double screen_y = graph_height * (-0.5*ndc_y + 0.5);
	return {(float)screen_x, (float)screen_y};
}

Coordinates2D* coordinates_2d(Function* func, double xmi, double xmx, Function* main=nullptr, double dx=0.001, double ymi=-INF, double ymx=INF)
{
	Coordinates2D* coords = new Coordinates2D();

	double range_min = -INF;
	double range_max = INF;

	std::vector<Vector2> xy_coords;
	std::vector<Vector2> screen_coords;
	for (double x = xmi; x <= xmx; x+=dx) {
		double y = ast_solve(func->function, x);
		if (y >= ymi && y <= ymx) {
			if (y<range_max) range_max = y;
			if (y>range_min) range_min = y;
			xy_coords.emplace_back(x, y);
		}
	}
	coords->xy_coordinates = xy_coords;
	coords->range = {range_max, range_min};
	coords->domain = !main ? std::pair<double, double>(xmi, xmx) : main->coords->domain;

	double domain = coords->domain.second - coords->domain.first;
	double range = coords->range.second - coords->range.first;
	double scale = 2.0 / std::max(domain, range);
	double x_offset = -(coords->domain.first + domain / 2.0) * scale;
	double y_offset = -(coords->range.first + range / 2.0) * scale;
	
	// AXES COORDS
	if (!main) {
		std::vector<Vector2> axes_coords = {
			{0, (float)coords->range.first},	
			{0, (float)coords->range.second},	
			{(float)coords->domain.first, 0},	
			{(float)coords->domain.second, 0}
		};
		for (const auto& point : axes_coords) {
			coords->sa_coords.emplace_back(normalize(scale, point.x, point.y, x_offset, y_offset));
		}
	}
	// FUNCTION COORDS
    for (const auto& point : xy_coords) {
        screen_coords.emplace_back(normalize(scale, point.x, point.y, x_offset, y_offset));
    }

	coords->screen_coordinates = screen_coords;
	func->coords = coords;

	// @DEBUG
	//For (screen_coords) {
	//	std::cout<<item.x<<" "<<item.y<<std::endl;
	//}
	return coords;
}

Function* create_function(String str)
{
	Function* fun = new Function();
	fun->str = str;
	fun->function = ast_tree(fun->str);
	return fun;
}

const std::vector<Color> colors = {
	{0,255,0,255},
	{255,0,0,255},
	{0,0,255,255},
	{255,0,255,255}
};

void draw_axes(Function* f)
{
	auto c = f->coords->sa_coords;
	DrawLineV(c[0], c[1], {214,214,214,255});
	DrawLineV(c[2], c[3], {214,214,214,255});
}

void draw(std::vector<Function*> f)
{
	for (int i = 0; i < f.size(); i++) {
		DrawLineStrip(f[i]->coords->screen_coordinates.data(), f[i]->coords->screen_coordinates.size(), colors[i]);
	}
}

std::vector<Function*> plot(std::vector<String> strs, double xmi, double ymi)
{
	std::vector<Function*> fs;
	Function* main = create_function(strs[0]);
	main->coords = coordinates_2d(main, xmi, ymi);
	fs.emplace_back(main);
	
	for (int i = 1; i < strs.size(); i++) {
		Function* f = create_function(strs[i]);
		f->coords = coordinates_2d(f, xmi, ymi, main);
		fs.emplace_back(f);
	}
	return fs;
}

// Allow negative numbers, variables, functions, etc.
// Variable before number (*)
// Draw Axes as a function?

int main()
{
	//auto fs = plot({"abs(x)^(2/3) + (9/10)*sin(20*abs(x)))*((3-(abs(x))^2)^(1/2))"}, -2, 5);
	auto fs = plot({"sinx"},-PI,PI);

	InitWindow(1400,1000, "Math.exe");
	while (!WindowShouldClose()) {
		ClearBackground(Color());
		BeginDrawing();
		//DrawLine(400, 0, 400, 800, {214,214,214,255});
		//DrawLine(0, 400, 800, 400, {214,214,214,255});
		draw(fs);
		draw_axes(fs[0]);
		EndDrawing();
	}
}
