#include "normal_form.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace nf;

int main() {
    // ===== 1. 构造 K_normal（只做一次） =====
    Poly H_complex = load_poly_from_file(
        "d:/keyan/projects/catalogue/Normal_Form/normal_form/H_complex_terms.txt");

    const Real lambda = 2.932055933522975L;
    const Real omegaY = 2.334385885011224L;
    const Real omegaZ = 2.268831094896147L;
    const std::array<Cplx, 3> nu = {
        Cplx{lambda, 0},
        Cplx{0,      omegaY},
        Cplx{0,      omegaZ}
    };

    const int N = 10; // 你生成 H_complex 用的阶数
    auto res = lie_transformation(H_complex, nu, N);
    res.K_normal.prune(1e-15L);

    // ===== 2. 打印复坐标到正规形坐标的转换公式 =====
    nf::Poly Q1, Q2, Q3, P1, P2, P3;
    build_coordinate_transforms_forward(res, N, Q1, Q2, Q3, P1, P2, P3);

    std::array<std::string, 6> names = { "q1","q2","q3","p1","p2","p3" };

    std::cout << "Q1(q,p) = " << Q1.to_string(names, 100) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "Q2(q,p) = " << Q2.to_string(names, 100) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "Q3(q,p) = " << Q3.to_string(names, 100) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "P1(q,p) = " << P1.to_string(names, 100) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "P2(q,p) = " << P2.to_string(names, 100) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "P3(q,p) = " << P3.to_string(names, 100) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    // ===== 3. 打印正规形表达式及作用角变量表达式 =====
    auto full = nf::extract_actions_expansion_full(res.K_normal, 6);

    std::cout << std::setprecision(9);
    std::cout << "K(I) = ";
    bool first = true;
    for (const auto& kv : full.terms) {
        const auto& eI = kv.first;
        const auto& c = kv.second;
        if (!first) std::cout << " + ";
        first = false;
        std::cout << "(" << c.real() << (c.imag() >= 0 ? "+" : "") << c.imag() << "i)*"
            << "I1^" << eI[0] << "*I2^" << eI[1] << "*I3^" << eI[2];
    }
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    Real tol = 1e-6L;

    std::cout << std::setprecision(9);
    std::cout << "K(I) = ";
    for (const auto& kv : full.terms) {
        const auto& eI = kv.first;
        const auto& c = kv.second;
        if (std::abs(c) < tol) continue;          // 小于阈值的项直接跳过

        if (!first) std::cout << " + ";
        first = false;

        Real re = c.real();
        Real im = c.imag();
        // 如果虚部很小，可以只打印实部
        if (std::abs(im) < tol) im = 0;

        std::cout << "(" << re << (im >= 0 ? "+" : "") << im << "i)*"
            << "I1^" << eI[0] << "*I2^" << eI[1] << "*I3^" << eI[2];
    }
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout << res.K_normal.to_string(
        { "q1","q2","q3","p1","p2","p3" }, 80) << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";


    // ===== 4. 读 MATLAB 导出的 qp_points.txt =====
    std::ifstream in("d:/keyan/projects/catalogue/qp_points.txt");
    if (!in) {
        std::cerr << "Cannot open qp_points.txt\n";
        return 1;
    }

    std::ofstream out("d:/keyan/projects/catalogue/action_angle.txt");
    out << std::setprecision(16);

    while (true) {
        long double qr1, qi1, qr2, qi2, qr3, qi3;
        long double pr1, pi1, pr2, pi2, pr3, pi3;
        if (!(in >> qr1 >> qi1 >> qr2 >> qi2 >> qr3 >> qi3
            >> pr1 >> pi1 >> pr2 >> pi2 >> pr3 >> pi3)) {
            break; // EOF
        }

        std::array<Cplx, 6> qp = {
            Cplx{qr1, qi1},
            Cplx{qr2, qi2},
            Cplx{qr3, qi3},
            Cplx{pr1, pi1},
            Cplx{pr2, pi2},
            Cplx{pr3, pi3}
        };

        ActionAngle aa = evaluate_action_angle(qp);

        // I1 I2 I3  phi1 phi2 phi3
        out << aa.I[0].real() << " " << aa.I[0].imag() << " "
            << aa.I[1].real() << " " << aa.I[1].imag() << " "
            << aa.I[2].real() << " " << aa.I[2].imag() << " "
            << aa.phi[0] << " " << aa.phi[1] << " " << aa.phi[2] << "\n";
    }

    return 0;
}

