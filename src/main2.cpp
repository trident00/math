
// Types and Data Structures
struct PlotData {
    std::string function;
    std::vector<Vector2> coordinates;
    std::vector<Vector2> screen_coordinates;
    bool smooth;
    bool incremental;
    double x_min, y_min, x_max, y_max;
    double x_scale, y_scale, x_offset, y_offset;
    int current_point = 0;
    bool drawing_complete = false;
};

// Evaluate Function and return STANDARD coordinates
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
void c_minmax(PlotData &data)
{
	data.x_min = std::numeric_limits<float>::max();
	data.x_max = std::numeric_limits<float>::lowest();
	data.y_min = std::numeric_limits<float>::max();
	data.y_max = std::numeric_limits<float>::lowest();

	for (const auto& vec : coordinates) {
		data.x_min = std::min((float)x_min, vec.x);
		data.x_max = std::max((float)x_max, vec.x);
		data.y_min = std::min((float)y_min, vec.y);
		data.y_max = std::max((float)y_max, vec.y);
	}
}
void c_scale(PlotData &data)
{
	double x_range = x_max - x_min;
	double y_range = y_max - y_min;
	// Normalization of coords
	double scale = std::max(x_range, y_range);
	x_scale = 2.0 / scale;
	y_scale = 2.0 / scale;
	x_offset = -(x_min + x_range / 2.0) * x_scale;
	y_offset = -(y_min + y_range / 2.0) * y_scale;
}
