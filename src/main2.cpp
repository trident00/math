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

double scale(std::pair<double,double> domain, double val)
{
	double x_center = (domain.second - domain.first) / 2;
	return val * (graph_width/2) / x_center + graph_width/2;
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
	For (xy_coords) {
		screen_coords.emplace_back(scale(coords->domain, item.x),scale(coords->domain, -item.y));
	}
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

void draw(Function* F, Color color)
{
	DrawLineStrip(F->coords->screen_coordinates.data(), F->coords->screen_coordinates.size(), color);
}

// Allow negative numbers, variables, functions, etc.
// Variable before number (*)
// Draw Axes as a function?

int main()
{
	auto* f = create_function("3xsinx");
	auto* f2 = create_function("-3x^2");
	auto* s = coordinates_2d(f, -10, 10);
	auto* s2 = coordinates_2d(f2, -2, 2);
	InitWindow(1000,800, "Math.exe");
	while (!WindowShouldClose()) {
		ClearBackground(Color());
		BeginDrawing();
		DrawLine(400, 0, 400, 800, {214,214,214,255});
		DrawLine(0, 400, 800, 400, {214,214,214,255});
		draw(f2, {255,0,0,255});
		EndDrawing();
	}
	//For (s) {
	//	std::cout<<item.x<<" "<<item.y<<std::endl;
	//}
}
