
#include <iostream>
#include <vector>
#include <functional>
#include <chrono>
#define _USE_MATH_DEFINES
#include <math.h>
#include <sstream>
using namespace std::chrono;

bool is_epsilon(const double& val, const double& compare_val, const double& eps)
{
    if (std::fabs(compare_val) <= std::numeric_limits<double>::epsilon())
        return std::fabs(val) < eps;

    return std::fabs(val / compare_val) < eps;
}

struct Scheme
{
    double m_ñoeff, m_tau, m_h;
    Scheme(double ñoeff, double tau, double h) : m_ñoeff(ñoeff), m_tau(tau), m_h(h) {}
};

struct CartesianScheme : Scheme
{
    double m_k;
    CartesianScheme(double ñoeff, double tau, double h) : Scheme(ñoeff, tau, h), m_k(ñoeff* tau / h / h) {}
    double operator()(double x, double y1, double y2, double y3)
    {
        return m_k * (y1 - 2 * y2 + y3) + y2;
    }
};

struct PolarScheme : Scheme
{
    double m_k1, m_k2;
    PolarScheme(double ñoeff, double tau, double h) : Scheme(ñoeff, tau, h), m_k1(ñoeff* tau / h / h), m_k2(ñoeff* tau / 2 / h) {}
    double operator()(double r, double y1, double y2, double y3)
    {
        return m_k1 * (y1 - 2 * y2 + y3) + m_k2 / r * (y3 - y1) + y2;
    }
};

// kU * dU/dx + kdU * U = C
struct BoundaryCondition3rd
{
    double m_kU, m_kdU, m_C;
    BoundaryCondition3rd(double kU, double kdU, double C) :m_kU(kU), m_kdU(kdU), m_C(C) {}

    double operator()(double h, double t_curr, double y1, double y2)
    {
        double divider = 2 * h * m_kU / m_kdU - 3;
        if (is_epsilon(divider, 3, 0.00001))
            return (m_C - y1) / (h * m_kU / m_kdU - 1); //O(h)

        return (2 * h / m_kdU * m_C - 4 * y1 + y2) / divider; //O(h*h)
    }
};

class ExplicitDifferenceSchemeException : public std::exception
{
public:
    ExplicitDifferenceSchemeException(char const* const msg) : std::exception(msg) {}
};

struct Timer
{
    decltype(std::chrono::steady_clock::now()) t;
    Timer()
    {
        t = std::chrono::steady_clock::now();
    }
    ~Timer()
    {
        std::cout << "Elapsed: " << static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t).count()) / 1000 / 1000 << "ms" << std::endl;
    }
};

double error(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size())
        throw std::exception("Vector sizes must be equal");

    double summ = 0;
    double compare = 0;
    for (decltype(v1.size()) i = 0; i < v1.size(); ++i)
    {
        double diff = (v1[i] - v2[i]);
        summ += diff * diff;
        compare += std::max(std::abs(v1[i]), std::abs(v2[i]));
    }

    if (compare <= std::numeric_limits<double>::epsilon())
        return summ;

    return summ / compare;
}

using bc_t = std::function<double(double, double, double, double)>;

template<class scheme_t>
std::vector<double> explicit_difference_scheme(const std::vector<double>& y0, double ñoeff, double a, double b, double t, bc_t bc_a, bc_t bc_b, double eps = 0.001, unsigned int max_iter = 100)
{
    if (y0.size() < 2)
        throw ExplicitDifferenceSchemeException("Minimum 2 points are required");

    if (!bc_a || !bc_b)
        throw ExplicitDifferenceSchemeException("Boundary conditions must be defined");

    if (b < a)
        throw ExplicitDifferenceSchemeException("Incorrect interval");

    if (ñoeff <= 0)
        throw ExplicitDifferenceSchemeException("Incorrect coefficient");

    if (t < 0)
        throw ExplicitDifferenceSchemeException("Incorrect duration");

    Timer timer;

    double h = (b - a) / (y0.size() - 1);
    if (h <= std::numeric_limits<double>::epsilon() * b)
        throw ExplicitDifferenceSchemeException("Interval is too small");

    unsigned int N = static_cast<unsigned int>(std::round(2 * ñoeff * t / h / h + 0.5));
    unsigned int iter = 1;
    double err = 0;
    std::remove_const<std::remove_reference<decltype(y0)>::type>::type y(y0.size(), 0), y_prev(y0.size(), 0);
    struct TmpVal
    {
        double val;
        TmpVal* next;
    } tmp[2] = { {0, tmp + 1}, {0, tmp} };

    do
    {
        y_prev = y;
        y = y0;
        double tau = t / N;
        std::cout << "iter=" << iter << ", h=" << h << ", tau=" << tau;
        scheme_t scheme(ñoeff, tau, h);
        for (decltype(N) n = 1; n <= N; ++n)
        {
            double t_curr = n * tau;

            //inner nodes
            TmpVal* curr_tmp = &tmp[0];
            if (y.size() > 2)
            {
                curr_tmp->val = scheme(a + h, y[0], y[1], y[2]);
                if (y.size() > 3)
                {
                    curr_tmp = &tmp[1];
                    curr_tmp->val = scheme(a + 2 * h, y[1], y[2], y[3]);
                }
            }

            for (decltype(y.size()) i = 3; i < y.size() - 1; ++i)
            {
                double x = a + h * i;
                curr_tmp = curr_tmp->next;
                y[i - 2] = curr_tmp->val;
                curr_tmp->val = scheme(x, y[i - 1], y[i], y[i + 1]);
            }
            if (y.size() > 2)
            {
                y[y.size() - 2] = curr_tmp->val;
                if (y.size() > 3)
                    y[y.size() - 3] = curr_tmp->next->val;
            }

            //boundary conditions
            y[0] = bc_a(h, t_curr, y[1], y[2]);
            y[y.size() - 1] = bc_b(h, t_curr, y[y.size() - 2], y[y.size() - 3]);
        }

        err = error(y, y_prev);
        std::cout << ", N = " << N << ", err = " << err << std::endl;
        N *= 2;

    } while (err > eps && iter++ <= max_iter);

    return y; //hope for RVO;
}

std::vector<double> test_ring(std::vector<double>& U, double R1, double R2, std::chrono::duration<double, std::ratio<1>> t, double density, double Cp, double thermal_conductivity, double heat_transfer_coefficient1, double heat_transfer_coefficient2, double T1, double T2, double eps = 0.1)
{
    double a = thermal_conductivity / Cp / density;
    BoundaryCondition3rd r1(-heat_transfer_coefficient1, 1, -T1 / heat_transfer_coefficient1);
    BoundaryCondition3rd r2(-heat_transfer_coefficient2, 1, -T2 / heat_transfer_coefficient2);
    return explicit_difference_scheme<PolarScheme>(U, a, R1, R2, t.count(), r1, r2, eps, 10);
}

void test1(int nodes, double eps);
void test2(int nodes, double eps);
void test3(int nodes, double eps);

std::ostream& operator<< (std::ostream& o, std::vector<double> v)
{
    for (auto a : v)
        o << a << ", ";

    return o;
}

int main(int argc, char* argv[])
{
    /* {
        Timer timer;
        std::vector<double> y(32, 10);
        double R1 = 1;
        double R2 = 10;
        auto t = 2s;
        double density = 1;
        double Cp = 1;
        double thermal_conductivity = 4;
        double heat_transfer_coefficient = 1;
        double T1 = 10;
        double T2 = 100;

        try
        {
            auto res = test_ring(y, R1, R2, t, density, Cp, thermal_conductivity, heat_transfer_coefficient, heat_transfer_coefficient, T1, T2, 0.001);
            std::cout << res << std::endl;
        }
        catch (const ExplicitDifferenceSchemeException& err)
        {
            std::cerr << err.what() << std::endl;
        }
        catch (const std::exception& err)
        {
            std::cerr << err.what() << std::endl;
        }
        catch (...)
        {
            std::cerr << "Unknown error" << std::endl;
        }
    }
    */

    std::string test = "test1";
    int nodes = 10;
    double eps = 0.001;

    if (argc > 1) {
        std::istringstream ss(argv[1]);
        ss >> test;
        if (argc > 2) {
            std::istringstream ss(argv[2]);
            ss >> nodes;
            if (argc > 3) {
                std::istringstream ss(argv[3]);
                ss >> eps;
            }
        }
    }

    if (test == "test1")
        test1(nodes, eps);
    else if (test == "test2")
        test2(nodes, eps);
    else if (test == "test3")
        test3(nodes, eps);
    else
        std::cout << "No such test" << std::endl;

    return 0;
}

//========================================================================================================================================

double integral(const std::function<double(double)>& f, double a, double b, double eps = 0.00001, int minN = 2, int maxN = 3000)
{
    double res = 0, res_prev = 0, summ = (f(a) + f(b)) / 2;

    int N = 1;
    do
    {
        double h = (b - a) / N;
        for (int i = 1; i <= N - 1; i += 2)
            summ += f(a + h * i);

        res_prev = res;
        res = summ * h;
        N *= 2;
    } while ((!is_epsilon(std::fabs(res - res_prev), std::fabs(res + res_prev) / 2, eps) || minN >= N) && N < maxN);


    return res;
}


// dU/dt = 4 * d2U/dx2 ; U(x,0)=1, U(0,t)=0, U(1,t)=1
double exact_solution_decart(double x, double t, double eps = 0.001)
{
    double summ = x;
    double add;
    int n = 1;
    do
    {
        //double tmp = M_PI * n;
        add = 2.0 / M_PI / n * exp(-4 * M_PI * M_PI * n * n * t) * sin(M_PI * n * x);
        summ += add;
        ++n;
    } while (std::abs(add) > eps);
    return summ;
}

// dU/dt = a d2U/dx2, U(0,t)=U(l,t)=0, u(x, 0)= T0(x)
double exact_solution_decart00(double x, double t, double a, double l, std::function<double(double)> T0, double eps = 0.0001)
{
    double summ = 0;
    double add;
    int n = 1;
    do
    {
        add = exp(-a * n * n * M_PI * M_PI * t / l / l) * sin(n * M_PI * x / l);
        double A = 2. / l * integral([n, l, T0](double xx) {return T0(xx) * sin(M_PI * n * xx / l); }, 0, l, eps);
        summ += add * A;
        ++n;
    } while (std::abs(add) > eps);

    return summ;
}


unsigned long double fact(int n)
{
    unsigned long double res = 1;
    for (int i = 1; i <= n; ++i)
        res *= i;

    return res;
};

template<int s>
double J(double x, double eps = 0.0001)
{
    double summ = 0, add;
    int n = 0;
    do
    {
        add = (n % 2 ? -1 : 1) * pow(x / 2, 2 * n + s) / fact(n) / fact(n + s);
        summ += add;
        ++n;
    } while (std::abs(add) > eps && n <= 100);
    return summ;
};

// dU/dt = 4 * (d2U/dr2 + 1/r dU/dr) ; 0 <= r < 8, 0 < t < T, U(r,0) = 64 - r^2, u(8,t) = 0;
double exact_solution_polar(double r, double t, double eps = 0.001)
{
    static std::vector<double> bessel = { 2.404, 5.520, 8.654, 11.792, 14.931, 18.076, 21.212, 24.353, 27.494 };
    auto tt = J<1>(18);
    double summ = 0;
    for (size_t i = 0; i < bessel.size(); ++i)
    {
        summ += 1. / bessel[i] / bessel[i] / bessel[i] / J<1>(bessel[i]) * exp(-bessel[i] * bessel[i] / 16 * t) * J<0>(bessel[i] / 8 * r);
    }
    summ *= 512;

    return summ;
}

void dump_result(const std::vector<double>& y, std::function<double(double)> exact_y)
{
    std::cout << "No\t<==>\texact\t<==>\tcalculated\tdiff=\t(%)" << std::endl;

    for (int i = 0; i < y.size(); ++i)
    {
        double x = 1. / (y.size() - 1) * i;
        //std::cout << x << ", ";
        double exact = exact_y(x);
        std::cout << i << "\t<==>\t" << exact << "\t<==>\t" << y[i] << ",\tdiff=" << exact - y[i] << "\t(" << std::abs(exact - y[i]) / std::max(std::abs(exact), std::abs(y[i])) * 100 << "%)" << std::endl;
    }
}

//TEST dekart  dU/dt = 4 * d2U/dt2 ; U(x,0)=1, U(0,t)=0, U(1,t)=1
void test1(int nodes, double eps)
{
    std::vector<double> y(nodes, 1.0);
    double t = 6;

    y = explicit_difference_scheme<CartesianScheme>(y, 4, 0, 1, t, [](double, double, double, double) {return 0; }, [](double, double, double, double) {return 1; }, eps);
    dump_result(y, [t](double x) {return exact_solution_decart(x, t);});
}

//TEST dekart  dU/dt = a * d2U/dt2 ; U(x,0)=T0(x), U(0,t)=U(l,t)=0
void test2(int nodes, double eps)
{
    std::vector<double> y(nodes);
    double l = 100;
    double a = 9;

    auto T0 = [l](double x)
    {
        return x < l / 2 ? 200 / l * x : -200 / l * x + 200;
        //return -(x - l / 2)*(x - l / 2) + 100;
    };

    for (int i = 0; i < y.size(); ++i)
    {
        double x = l / (y.size() - 1) * i;
        y[i] = T0(x);
    }

    double t = 1;
    y = explicit_difference_scheme<CartesianScheme>(y, a, 0, l, t, [](double, double, double, double) {return 0; }, [](double, double, double, double) {return 0; }, eps);

    dump_result(y, [t, a, l, T0](double x) {return exact_solution_decart00(l * x, t, a, l, T0); });
}

//TEST polar dU/dt = 4 * (d2U/dr2 + 1/r dU/dr) ; 0 <= r < 8, 0 < t < T, U(r,0) = 64 - r^2, u(8,t) = 0;
void test3(int nodes, double eps)
{
    std::vector<double> y(nodes, 0.0);
    for (int i = 0; i < y.size(); ++i)
    {
        double x = 8. / (y.size() - 1) * i;
        y[i] = 64 - x * x;
    }

    double t = 5;
    auto bc_a = [](double h, double t_curr, double y1, double y2)
    {
        return (4 * y1 - y2) / 3;
    };
    y = explicit_difference_scheme<PolarScheme>(y, 4, 0, 8, t, bc_a, [](double, double, double, double) {return 0; }, eps);

    dump_result(y, [t](double x) {return exact_solution_polar(8 * x, t); });
}
