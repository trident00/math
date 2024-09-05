#include <iostream>
#include <chrono>
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


