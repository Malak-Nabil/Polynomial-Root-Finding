#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>


using namespace std;


int validity_order_check(string equation) {
	//Checks order of the equation
	int p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, f = 0;
	for (int i = 0; i < equation.size(); i++) {
		if (equation[i] == 'x' && (equation[i + 1] != '^' || i == equation.size() - 1)) {
			p1 = i + 1;
		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '2') {
			p2 = i + 1;
		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '3') {
			p3 = i + 1;
		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '4') {
			p4 = i + 1;
		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '5') {
			p5 = i + 1;
		}

	}

	if (p1 > 0 && p2 == 0 && p3 == 0 && p4 == 0 && p5 == 0) {
		cout << "The equation is from the first degree" << endl;
		f = 1;
	}

	if (p2 > 0 && p3 == 0 && p4 == 0 && p5 == 0) {
		cout << "The equation is from the second degree" << endl;
		f = 2;

		if (p1 > p2) {
			cout << "Equation is in right arrangement" << endl;
		}
		else {
			cout << "The arrangement is wrong, please enter it again" << endl;
		}
	}
	if (p3 > 0 && p4 == 0 && p5 == 0) {
		cout << "The equation is from the third degree" << endl;
		f = 3;

		if (p1 > p2 && p2 > p3) {
			cout << "Equation is in right arrangement" << endl;
		}
		else {
			cout << "The arrangement is wrong, please enter it again" << endl;
		}
	}
	if (p4 > 0 && p5 == 0) {
		cout << "The equation is from the fourth degree" << endl;
		f = 4;
		if (p1 > p2 && p2 > p3 && p3 > p4) {

			cout << "Equation is in right arrangement" << endl;
		}
		else {
			cout << "The arrangement is wrong, please enter it again" << endl;
		}
	}
	if (p5 > 0) {
		cout << "The equation is from the fifth degree" << endl;
		f = 5;
		if (p1 > p2 && p2 > p3 && p3 > p4 && p4 > p5) {
			cout << "Equation is in right arrangement" << endl;
		}
		else {
			cout << "The arrangement is wrong, please enter it again" << endl;
		}
	}
	//Checks that equation format is correct
	for (int i = 0; i < equation.size(); i++) {
		//Checks that a number exists after + or - sign
		if ((equation[i] == '+' || equation[i] == '-') && (equation[i + 2] == ' '
			|| equation[equation.size() - 1] == '+' || equation[equation.size() - 1] == '-')) {
			cout << "The form is wrong please enter it again" << endl;
			f = 0;
			break;
		}
		//Checks that power of x exists and is a number between 1 and 5
		else if (equation[i] == '^' && !(equation[i + 1] > 49 && equation[i + 1] < 54)) {
			cout << "The form is wrong please enter it again" << endl;
			f = 0;
			break;
		}
		//Checks that power sign exists
		else if (equation[i] == 'x' && equation[i + 1] != '^' && equation[i + 1] != '+' && equation[i + 1] != '-') {
			cout << "The form is wrong please enter it again" << endl;
			f = 0;
			break;
		}
	}
	return f;
}

//Solving first degree equation
void first_degree(double& x_0, double& x_1, string equation, vector<double>& coeff) {
	string s_0, s_1;
	int f0 = 0, f1 = 0;
	double a, b, solution1;
	for (int i = 0; i < equation.size(); i++) {
		if (equation[i] == 'x' && equation[i + 1] != '^') {
			f1 = i;
			s_1 = equation.substr(0, f1); //gets coefficient of x^1
			x_1 = stoi(s_1); //converts string to int
			coeff.push_back(x_1); //stores coefficient of x^1 in vector
			f0 = i + 1;
			s_0 = equation.substr(f0, equation.size() - f0); //gets coefficient of x^0
			x_0 = stoi(s_0); //converts string to int
			coeff.push_back(x_0); //stores coefficient of x^0 in vector
		}

	}
	a = coeff[0]; //coefficient of x^1
	b = coeff[1]; //coefficient of x^0
	solution1 = -b / a;
	cout << "Solution = " << solution1;

}

//Solving second degree equation
void second_degree(double& x_0, double& x_1, double& x_2, string equation, vector<double>& coeff) {
	string s_0, s_1, s_2;
	int f0 = 0, f1 = 0, f2 = 0;
	double a, b, c, same, imaginary, delta;
	double solution1, solution2;


	for (int i = 0; i < equation.size(); i++) {
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '2') {
			f2 = i;
			s_2 = equation.substr(0, f2); //gets coefficient of x^2
			x_2 = stoi(s_2); //converts string to int
			coeff.push_back(x_2); //stores coefficient of x^2 in vector
		}
		if (equation[i] == 'x' && equation[i + 1] != '^') {
			f1 = i;
			s_1 = equation.substr(f2 + 3, f1 - f2 - 3); //gets coefficient of x^1
			x_1 = stoi(s_1); //converts string to int
			coeff.push_back(x_1); //stores coefficient of x^1 in vector
			f0 = i + 1;
			s_0 = equation.substr(f0, equation.size() - f0); //gets coefficient of x^0
			x_0 = stoi(s_0); //converts string to int
			coeff.push_back(x_0); //stores coefficient of x^0 in vector
		}

	}
	a = coeff[0]; //coefficient of x^2
	b = coeff[1]; //coefficient of x^1
	c = coeff[2]; //coefficient of x^0

	delta = pow(b, 2) - 4 * a * c;
	same = (-b) / (2 * a);
	solution1 = (-(b)+(sqrt(delta)) / (2 * a));
	solution2 = (-(b)-(sqrt(delta)) / (2 * a));

	cout << "Solution is presented as the following shape (Real No , Imaginary No)" << endl;
	//If delta = 0, then the 2 solutions are equal and real
	if (delta == 0)
	{
		cout << "solution1 = solution2 = " << (complex<double>)same << endl;
	}
	//If delta is negative,there are imaginary roots
	else if (delta < 0)
	{
		imaginary = sqrt(-delta) / (2 * a);
		cout << "solution 1 = " << same << "+" << imaginary << "i" << endl;
		cout << "solution 2 = " << same << "-" << imaginary << "i" << endl;
	}
	//If delta is positive, then the 2 solutions are real
	else
	{
		cout << "solution1=" << (complex<double>)solution1 << endl;
		cout << "solution2=" << (complex<double>)solution2 << endl;
	}

}

//Solving third degree equation
void third_degree(double& x_0, double& x_1, double& x_2, double& x_3, string equation, vector<double>& coeff) {
	string s_0, s_1, s_2, s_3;
	int f0 = 0, f1 = 0, f2 = 0, f3 = 0;
	double a, b, c, d;
	const double PI = 3.14159265359;
	double solution1, solution2, solution3;
	for (int i = 0; i < equation.size(); i++) {
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '3') {
			f3 = i;
			s_3 = equation.substr(0, f3); //gets coefficient of x^3
			x_3 = stoi(s_3); //converts string to int
			coeff.push_back(x_3); //stores coefficient of x^3 in vector
		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '2') {
			f2 = i;
			s_2 = equation.substr(f3 + 3, f2 - f3 - 3); //gets coefficient of x^2
			x_2 = stoi(s_2); //converts string to int
			coeff.push_back(x_2); //stores coefficient of x^2 in vector
		}
		if (equation[i] == 'x' && equation[i + 1] != '^') {
			f1 = i;
			s_1 = equation.substr(f2 + 3, f1 - f2 - 3); //gets coefficient of x^1
			x_1 = stoi(s_1); //converts string to int
			coeff.push_back(x_1); //stores coefficient of x^1 in vector
			f0 = i + 1;
			s_0 = equation.substr(f0, equation.size() - f0); //gets coefficient of x^0
			x_0 = stoi(s_0); //converts string to int
			coeff.push_back(x_0); //stores coefficient of x^0 in vector
		}

	}
	a = coeff[0]; //coefficient of x^3
	b = coeff[1]; //coefficient of x^2
	c = coeff[2]; //coefficient of x^1
	d = coeff[3]; //coefficient of x^0
	double p = (pow(b, 2) - (3.0 * a * c)) / (9.0 * pow(a, 2));
	double q = (9.0 * a * b * c - 27.0 * pow(a, 2) * d - 2.0 * pow(b, 3)) / (54.0 * pow(a, 3));
	double Z = b / (3.0 * a);
	double disc = (pow(p, 3)) - (pow(q, 2));

	//If disc is positive, the three solutions are real
	if (disc > 0)
	{
		double alpha = acos(q / (p * sqrt(p)));
		double r = 2.0 * sqrt(p);
		solution1 = r * cos((alpha + 2.0 * 0 * PI) / 3.0) - Z;
		solution2 = r * cos((alpha + 2.0 * 1 * PI) / 3.0) - Z;
		solution3 = r * cos((alpha + 2.0 * 2 * PI) / 3.0) - Z;
		cout << "solution 1 = " << solution1 << endl;
		cout << "solution 2 = " << solution2 << endl;
		cout << "solution 3 = " << solution3 << endl;
	}
	else
	{
		double Beta = cbrt(q + sqrt(-disc));
		double gamma = cbrt(q - sqrt(-disc));
		solution1 = Beta + gamma - Z;
		cout << "Solution 1 = " << solution1 << endl;


		double RealNo = -0.5 * (Beta + gamma) - Z;
		double ImagNo = (Beta - gamma) * sqrt(3.0) / 2.0;

		//If disc=0, then the three solutions are real and two of the 3 solutions are equal
		if (disc == 0.0)
		{
			cout << "solution 2 = " << RealNo << endl;
			cout << "solution 3 = " << RealNo << endl;
		}

		//If disc i negative, the first solution is real and the second and third solutions are imaginary
		else
		{
			cout << "solution 2 = " << RealNo << " + " << ImagNo << " i " << endl;
			cout << "solution 3 = " << RealNo << " - " << ImagNo << " i " << endl;
		}
	}

}

//Solving fourth degree equation
void fourth_degree(double& x_0, double& x_1, double& x_2, double& x_3, double& x_4, string equation, vector<double>& coeff) {
	string s_0, s_1, s_2, s_3, s_4;
	int f0 = 0, f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	double a, b, c, d, e;
	complex<double> p1, p2, p3, p4, p5, p6, solution1, solution2, solution3, solution4;
	for (int i = 0; i < equation.size(); i++) {
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '4') {
			f4 = i;
			s_4 = equation.substr(0, f4); //gets coefficient of x^4
			x_4 = stoi(s_4); //converts string to int
			coeff.push_back(x_4); //stores coefficient of x^4 in vector

		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '3') {
			f3 = i;
			s_3 = equation.substr(f4 + 3, f3 - f4 - 3); //gets coefficient of x^3
			x_3 = stoi(s_3); //converts string to int
			coeff.push_back(x_3); //stores coefficient of x^3 in vector
		}
		if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '2') {
			f2 = i;
			s_2 = equation.substr(f3 + 3, f2 - f3 - 3); //gets coefficient of x^2
			x_2 = stoi(s_2); //converts string to int
			coeff.push_back(x_2); //stores coefficient of x^2 in vector
		}
		if (equation[i] == 'x' && equation[i + 1] != '^') {
			f1 = i;
			s_1 = equation.substr(f2 + 3, f1 - f2 - 3); //gets coefficient of x^1
			x_1 = stoi(s_1); //converts string to int
			coeff.push_back(x_1); //stores coefficient of x^1 in vector
			f0 = i + 1;
			s_0 = equation.substr(f0, equation.size() - f0); //gets coefficient of x^0
			x_0 = stoi(s_0); //converts string to int
			coeff.push_back(x_0); //stores coefficient of x^0 in vector
		}
	}
	a = coeff[0]; //coefficient of x^4
	b = coeff[1]; //coefficient of x^3
	c = coeff[2]; //coefficient of x^2
	d = coeff[3]; //coefficient of x^1
	e = coeff[4]; //coefficient of x^0

	p1 = 2.0 * c * c * c - 9.0 * b * c * d + 27.0 * a * d * d + 27.0 * b * b * e - 72.0 * a * c * e;
	p2 = p1 + sqrt(complex<double>(-4.0 * (pow(c * c - 3.0 * b * d + 12.0 * a * e, 3.0)) + p1 * p1));
	p3 = (pow(c, 2.0) - 3.0 * b * d + 12.0 * a * e) / (3.0 * a * pow((p2 / 2.0), (1.0 / 3.0))) + (pow((p2 / 2.0), (1.0 / 3.0))) / (3.0 * a);
	p4 = sqrt(complex<double>((pow(b, 2.0)) / (4.0 * pow(a, 2)) - ((2.0 * c) / (3.0 * a)) + p3));
	p5 = pow(b, 2.0) / (2.0 * pow(a, 2.0)) - (4.0 * c) / (3.0 * a) - p3;
	p6 = (-(pow(b, 3.0) / pow(a, 3.0)) + (4.0 * b * c) / pow(a, 2.0) - ((8.0 * d) / a)) / (4.0 * p4);

	solution1 = -(b / (4.0 * a)) - (p4 / 2.0) - (sqrt(complex<double>(p5 - p6))) / 2.0;
	solution2 = -(b / (4.0 * a)) - (p4 / 2.0) + (sqrt(complex<double>(p5 - p6))) / 2.0;
	solution3 = -(b / (4.0 * a)) + (p4 / 2.0) - (sqrt(complex<double>(p5 + p6))) / 2.0;
	solution4 = -(b / (4.0 * a)) + (p4 / 2.0) + (sqrt(complex<double>(p5 + p6))) / 2.0;

	//There are 4 solutions in form of real no + or - the imaginary no
	cout << "The roots are in shape of (real_part, imaginary_part):" << endl;
	cout << "Solution 1 = " << solution1 << endl
		<< "Solution 2 = " << solution2 << endl
		<< "Solution 3 = " << solution3 << endl
		<< "Solution 4 = " << solution4 << endl;
}

//Solving fifth degree equation
void fifth_degree(double& x_0, double& x_1, double& x_2, double& x_3, double& x_4, double& x_5, string equation, vector<double>& coeff) {
	int f0 = 0, f1 = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0;
	string s_0, s_1, s_2, s_3, s_4, s_5;
	vector<double>A;//first term
	vector<double>B;// last term
	vector<double>number;
	complex<double> p1, p2, p3, p4, p5, p6, solution2=0, solution3=0, solution4=0, solution5=0;
	double solution1=0;
	int a, b, c, d, e, f;

	for (int i = 0; i < equation.size(); i++) {
		if (equation[i] == 'x' && equation[i+1]=='^' && equation[i + 2] == '5') {
			f5 = i;
			s_5 = equation.substr(0, f5);
			x_5 = stoi(s_5);
			coeff.push_back(x_5);

		}
		else if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '4') {
			f4 = i;
			s_4 = equation.substr(f5 + 3, f4 - f5 - 3);
			x_4 = stoi(s_4);
			coeff.push_back(x_4);

		}
		else if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '3') {
			f3 = i;
			s_3 = equation.substr(f4 + 3, f3 - f4 - 3);
			x_3 = stoi(s_3);
			coeff.push_back(x_3);
		}
		else if (equation[i] == 'x' && equation[i + 1] == '^' && equation[i + 2] == '2') {
			f2 = i;
			s_2 = equation.substr(f3 + 3, f2 - f3 - 3);
			x_2 = stoi(s_2);
			coeff.push_back(x_2);
		}
		else if (equation[i] == 'x' && equation[i + 1] != '^') {
			f1 = i;
			s_1 = equation.substr(f2 + 3, f1 - f2 - 3);
			x_1 = stoi(s_1);
			coeff.push_back(x_1);
			f0 = i + 1;
			s_0 = equation.substr(f0, equation.size() - f0);
			x_0 = stoi(s_0);
			coeff.push_back(x_0);

		}
	}

	a = coeff[0]; //coefficient of x^5
	b = coeff[1]; //coefficient of x^4
	c = coeff[2]; //coefficient of x^3
	d = coeff[3]; //coefficient of x^2
	e = coeff[4]; //coefficient of x^1
	f = coeff[5]; //coefficient of x^0
	vector<int> fifth = { a,b,c,d,e,f };
	for (int y = 1; y <= a; y++)
	{
		if (a % y == 0)
		{
			A.push_back(y);
		}
	}
	for (int w = 1; w <= f; w++)
	{
		if (f % w == 0)
		{
			B.push_back(w);
		}
	}
	for (int i = 0; i < B.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
		{
			double z = B[i] / A[j];
			number.push_back(z);
			number.push_back(-z);
		}
	}
	for (int i = 0; i < number.size(); i++)
	{
		solution1 = number[i];
		double m = (a * pow(solution1, 5)) + (b * pow(solution1, 4)) + (c * pow(solution1, 3)) + (d * pow(solution1, 2)) + (e * solution1) + f;
		if (m == 0)
		{
			solution1 = number[i];
			break;

		}
	}
	vector <int> new_coef(5);
	new_coef[0] = solution1 * fifth[0];
	for (int i = 0; i < 4; i++) {
		new_coef[i + 1] = new_coef[i] * solution1 + fifth[i + 1];
	}
	a = new_coef[0];
	b = new_coef[1];
	c = new_coef[2];
	d = new_coef[3];
	e = new_coef[4];

	p1 = 2.0 * c * c * c - 9.0 * b * c * d + 27.0 * a * d * d + 27.0 * b * b * e - 72.0 * a * c * e;
	p2 = p1 + sqrt(complex<double>(-4.0 * (pow(c * c - 3.0 * b * d + 12.0 * a * e, 3.0)) + p1 * p1));
	p3 = (pow(c, 2.0) - 3.0 * b * d + 12.0 * a * e) / (3.0 * a * pow((p2 / 2.0), (1.0 / 3.0))) + (pow((p2 / 2.0), (1.0 / 3.0))) / (3.0 * a);
	p4 = sqrt(complex<double>((pow(b, 2.0)) / (4.0 * pow(a, 2)) - ((2.0 * c) / (3.0 * a)) + p3));
	p5 = pow(b, 2.0) / (2.0 * pow(a, 2.0)) - (4.0 * c) / (3.0 * a) - p3;
	p6 = (-(pow(b, 3.0) / pow(a, 3.0)) + (4.0 * b * c) / pow(a, 2.0) - ((8.0 * d) / a)) / (4.0 * p4);

	solution2 = -(b / (4.0 * a)) - (p4 / 2.0) - (sqrt(complex<double>(p5 - p6))) / 2.0;
	solution3 = -(b / (4.0 * a)) - (p4 / 2.0) + (sqrt(complex<double>(p5 - p6))) / 2.0;
	solution4 = -(b / (4.0 * a)) + (p4 / 2.0) - (sqrt(complex<double>(p5 + p6))) / 2.0;
	solution5 = -(b / (4.0 * a)) + (p4 / 2.0) + (sqrt(complex<double>(p5 + p6))) / 2.0;

	cout << "The roots are in shape of (real_part, imaginary_part):" << endl;
	cout << "Solution 1 = " << (complex<double>)solution1 << endl
		<< "Solution 2 = " << solution2 << endl
		<< "Solution 3 = " << solution3 << endl
		<< "Solution 4 = " << solution4 << endl
		<<"Solution 5 = " << solution5 << endl;
}



int main() {
	string equation;
	cout << "Please enter your equation in the form of a(n)x^(n)+a(n-1)x^(n-1)+...+a1x+a0 = 0 : " << endl;
	cout << "if there is a term non existing put its coeffecient =0 (0x^n)" << endl;
	cin >> equation; //user inputs the equation
	int degree;
	vector<double> coeff;
	double x_0 = 0, x_1 = 0, x_2 = 0, x_3 = 0, x_4 = 0, x_5 = 0;

	degree = validity_order_check(equation); //calling function to check format and order of equation

	//Determines the degree of the equation and calls its function to solve it
	if (degree == 5)
	{
		fifth_degree(x_0, x_1, x_2, x_3, x_4, x_5, equation, coeff);
	}
	if (degree == 4) {
		fourth_degree(x_0, x_1, x_2, x_3, x_4, equation, coeff);
	}
	if (degree == 3) {
		third_degree(x_0, x_1, x_2, x_3, equation, coeff);
	}
	if (degree == 2) {
		second_degree(x_0, x_1, x_2, equation, coeff);
	}
	if (degree == 1) {
		first_degree(x_0, x_1, equation, coeff);
	}

}