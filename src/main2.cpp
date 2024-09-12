#include "parser.h"
#include "raylib.h"
#include "ast.h"
#include <iostream>
#include <stdio.h>
#include <limits>

#define E		2.7182818284590452353602874713527
#define INF     std::numeric_limits<double>::max()

const int graph_width = 800;
const int graph_height = 800;

struct Coordinates2D
{
	std::vector<Vector2> xy_coordinates;
	std::pair<double, double> domain;
	std::pair<double, double> range;
	std::vector<Vector2> screen_coordinates;
	std::vector<Vector2> axes;

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

//double scale(std::pair<double,double> domain, double val)
//{
//	double x_center = (domain.second - domain.first) / 2;
//	return val * (graph_width/2) / x_center + graph_width/2;
//}
	
// RANGE>DOMAIN FUNCTIONS BROKEN
double scale(std::pair<double,double> domain, std::pair<double,double> range, double val)
{
	//double x_center = (domain.second - domain.first) / 2;
	//double y_center = (range.second - range.first) / 2;
	//double scale = std::max(x_center, y_center);
	double scale = range.first + ((val - domain.first) / (domain.second - domain.first)) * (range.second - range.first);

	return scale;
}

Coordinates2D* coordinates_2d(Function* func, double xmi, double xmx, Function* main=nullptr, double dx=0.01, double ymi=-INF, double ymx=INF)
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
	if (coords->domain > coords->range) {
		
	}
// Now, normalize and transform the coordinates to screen space
    for (const auto& point : xy_coords) {
        // Normalize x and y
        double x_norm = (point.x - xmi) / (coords->domain.second - coords->domain.first);
        double y_norm = (point.y - range_min) / (coords->range.second - coords->range.first);

        // Transform to screen space
        double screen_x = x_norm * x_scale * graph_width;
        double screen_y = (1 - y_norm) * y_scale * graph_height;  // Inverted y axis for screen

        screen_coords.emplace_back(screen_x, screen_y);
    }

	//For (xy_coords) {
	//	screen_coords.emplace_back(scale(coords->domain, coords->range, item.x),scale(coords->domain, coords->range, -item.y));
	//}
	coords->screen_coordinates = screen_coords;
	func->coords = coords;
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

void draw(std::vector<Function*> f)
{
	int i = 0;
	For (f) {
		DrawLineStrip(item->coords->screen_coordinates.data(), item->coords->screen_coordinates.size(), colors[i]);
		i++;
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

	InitWindow(1000,800, "Math.exe");
	while (!WindowShouldClose()) {
		ClearBackground(Color());
		BeginDrawing();
		DrawLine(400, 0, 400, 800, {214,214,214,255});
		DrawLine(0, 400, 800, 400, {214,214,214,255});
		draw(fs);
		EndDrawing();
	}
}
