#include "normal_form.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <fstream>

namespace nf {

    std::size_t Exp6Hash::operator()(const Exp6& e) const noexcept {
        // FNV-1a style hash
        std::size_t h = 1469598103934665603ull;
        for (std::uint8_t v : e) {
            h ^= static_cast<std::size_t>(v);
            h *= 1099511628211ull;
        }
        return h;
    }

    Poly::Poly(Cplx c) { add_term(Exp6{ 0, 0, 0, 0, 0, 0 }, c); }

    void Poly::add_term(const Exp6& e, Cplx c) {
        if (c == Cplx{ 0, 0 }) return;
        auto it = terms_.find(e);
        if (it == terms_.end()) {
            terms_.emplace(e, c);
        }
        else {
            it->second += c;
            if (it->second == Cplx{ 0, 0 }) terms_.erase(it);
        }
    }

    void Poly::add_poly(const Poly& other, Cplx scale) {
        if (scale == Cplx{ 0, 0 }) return;
        for (const auto& kv : other.terms_) add_term(kv.first, scale * kv.second);
    }

    Poly Poly::scaled(Cplx s) const {
        Poly out;
        if (s == Cplx{ 0, 0 }) return out;
        out.terms_.reserve(terms_.size());
        for (const auto& kv : terms_) out.terms_.emplace(kv.first, s * kv.second);
        return out;
    }

    Poly Poly::truncated(int max_degree) const {
        Poly out;
        out.terms_.reserve(terms_.size());
        for (const auto& kv : terms_) {
            if (total_degree(kv.first) <= max_degree) out.terms_.emplace(kv.first, kv.second);
        }
        return out;
    }

    Poly Poly::homogeneous(int degree) const {
        Poly out;
        for (const auto& kv : terms_) {
            if (total_degree(kv.first) == degree) out.terms_.emplace(kv.first, kv.second);
        }
        return out;
    }

    void Poly::prune(Real abs_threshold) {
        if (abs_threshold <= 0) return;
        for (auto it = terms_.begin(); it != terms_.end();) {
            if (std::abs(it->second) < abs_threshold) {
                it = terms_.erase(it);
            }
            else {
                ++it;
            }
        }
    }
    Poly load_poly_from_file(const std::string& path) {
        Poly H;
        std::ifstream in(path);
        if (!in) {
            throw std::runtime_error("Cannot open file: " + path);
        }
        int e1, e2, e3, e4, e5, e6;
        long double re, im;
        while (in >> e1 >> e2 >> e3 >> e4 >> e5 >> e6 >> re >> im) {
            Exp6 e = {
                static_cast<std::uint8_t>(e1),
                static_cast<std::uint8_t>(e2),
                static_cast<std::uint8_t>(e3),
                static_cast<std::uint8_t>(e4),
                static_cast<std::uint8_t>(e5),
                static_cast<std::uint8_t>(e6)
            };
            Cplx c{ re, im };
            H.add_term(e, c);
        }
        return H;
    }

    static std::string cplx_to_string(const Cplx& c) {
        std::ostringstream os;
        os << std::setprecision(18);
        const auto re = c.real();
        const auto im = c.imag();
        if (im == 0) {
            os << re;
        }
        else if (re == 0) {
            os << im << "i";
        }
        else {
            os << "(" << re << (im >= 0 ? "+" : "") << im << "i)";
        }
        return os.str();
    }

    std::string Poly::to_string(const std::array<std::string, 6>& var_names,
        int max_terms) const {
        std::vector<std::pair<Exp6, Cplx>> items;
        items.reserve(terms_.size());
        for (const auto& kv : terms_) items.push_back(kv);

        std::sort(items.begin(), items.end(),
            [](const auto& a, const auto& b) {
                const int da = total_degree(a.first);
                const int db = total_degree(b.first);
                if (da != db) return da > db;
                return a.first < b.first;
            });

        std::ostringstream os;
        int shown = 0;
        for (const auto& kv : items) {
            if (max_terms > 0 && shown >= max_terms) break;
            const auto& e = kv.first;
            const auto& c = kv.second;
            if (shown > 0) os << " + ";
            os << cplx_to_string(c);
            for (int i = 0; i < 6; i++) {
                const auto pow = e[i];
                if (pow == 0) continue;
                os << "*" << var_names[i];
                if (pow != 1) os << "^" << int(pow);
            }
            shown++;
        }
        if (terms_.empty()) return "0";
        if (max_terms > 0 && int(items.size()) > max_terms) {
            os << " + ...(" << (items.size() - max_terms) << " more terms)";
        }
        return os.str();
    }

    // Derivative of a monomial term: c * x^e
    // with respect to var idx -> (c*e[idx]) * x^(e - 1 at idx) if e[idx]>0
    static bool deriv_term(const Exp6& e, const Cplx& c, int idx, Exp6& out_e, Cplx& out_c) {
        const auto pow = e[idx];
        if (pow == 0) return false;
        out_e = e;
        out_e[idx] = static_cast<std::uint8_t>(pow - 1);
        out_c = c * Cplx{ static_cast<Real>(pow), 0 };
        return true;
    }

    static Exp6 add_exps(const Exp6& a, const Exp6& b) {
        Exp6 out = a;
        for (int i = 0; i < 6; i++) out[i] = static_cast<std::uint8_t>(a[i] + b[i]);
        return out;
    }

    Poly poisson_bracket(const Poly& f, const Poly& g) {
        // {f,g} = sum_j df/dqj * dg/dpj - df/dpj * dg/dqj
        Poly out;

        // Pre-split terms for g derivatives to reduce repeated work per f term is possible,
        // but we keep it simple first.
        for (int j = 0; j < 3; j++) {
            const int q_idx = j;       // 0..2
            const int p_idx = 3 + j;   // 3..5

            // Sum over monomials in f and g: build df/dq * dg/dp
            for (const auto& ft : f.terms()) {
                Exp6 ef1;
                Cplx cf1;
                if (!deriv_term(ft.first, ft.second, q_idx, ef1, cf1)) continue;
                for (const auto& gt : g.terms()) {
                    Exp6 eg1;
                    Cplx cg1;
                    if (!deriv_term(gt.first, gt.second, p_idx, eg1, cg1)) continue;
                    out.add_term(add_exps(ef1, eg1), cf1 * cg1);
                }
            }

            // subtract df/dp * dg/dq
            for (const auto& ft : f.terms()) {
                Exp6 ef2;
                Cplx cf2;
                if (!deriv_term(ft.first, ft.second, p_idx, ef2, cf2)) continue;
                for (const auto& gt : g.terms()) {
                    Exp6 eg2;
                    Cplx cg2;
                    if (!deriv_term(gt.first, gt.second, q_idx, eg2, cg2)) continue;
                    out.add_term(add_exps(ef2, eg2), -cf2 * cg2);
                }
            }
        }

        return out;
    }

    Poly compute_generator(const Poly& Hn, const std::array<Cplx, 3>& nu, Real denom_tol) {
        Poly Gn;
        for (const auto& kv : Hn.terms()) {
            const auto& e = kv.first;
            const auto& coeff = kv.second;

            const Exp6& k = e;
            const bool resonant_simple =
                (k[0] == k[3]) && (k[1] == k[4]) && (k[2] == k[5]);
            if (resonant_simple) continue;

            const int dq1 = int(k[3]) - int(k[0]);
            const int dq2 = int(k[4]) - int(k[1]);
            const int dq3 = int(k[5]) - int(k[2]);
            const Cplx denom = Cplx{ Real(dq1), 0 } *nu[0] + Cplx{ Real(dq2), 0 } *nu[1] +
                Cplx{ Real(dq3), 0 } *nu[2];
            if (std::abs(denom) <= denom_tol) continue;

            const Cplx gn_coeff = -coeff / denom;
            Gn.add_term(e, gn_coeff);
        }

        return Gn;
    }

    Poly extract_resonant_terms(const Poly& K_current, const std::array<Cplx, 3>& nu,
        int order_n, Real resonance_tol) {
        Poly Hn = K_current.homogeneous(order_n);
        Poly out;

        for (const auto& kv : Hn.terms()) {
            const auto& k = kv.first;
            const auto& coeff = kv.second;

            const int dq1 = int(k[3]) - int(k[0]);
            const int dq2 = int(k[4]) - int(k[1]);
            const int dq3 = int(k[5]) - int(k[2]);
            const Cplx inner = Cplx{ Real(dq1), 0 } *nu[0] + Cplx{ Real(dq2), 0 } *nu[1] +
                Cplx{ Real(dq3), 0 } *nu[2];
            if (std::abs(inner) <= resonance_tol) out.add_term(k, coeff);
        }

        return out;
    }

    Poly apply_lie_transform_H(const Poly& K, const Poly& G, int max_order, int n) {
        Poly out = K;
        Poly Lprev = poisson_bracket(K, G).truncated(max_order); // L_G K
        int k = 1;
        Real factorial = 1.0L;

        while (!Lprev.empty()) {
            // 估算这一阶的最小次数：2 + k*(n-2)
            int min_deg = 2 + k * (n - 2);
            if (min_deg > max_order) break;

            factorial *= static_cast<Real>(k);
            out.add_poly(Lprev, Cplx{ 1.0L / factorial, 0 });

            // 下一次括号
            Poly Lnext = poisson_bracket(Lprev, G).truncated(max_order);
            Lprev = std::move(Lnext);
            ++k;
        }

        return out.truncated(max_order);
    }

    NormalFormResult lie_transformation(const Poly& H_complex, const std::array<Cplx, 3>& nu,
        int max_order, Real denom_tol, Real resonance_tol) {
        NormalFormResult res;
        res.G.resize(max_order + 1);

        Poly K_current = H_complex;
        res.K_normal = K_current.homogeneous(2);

        for (int n = 3; n <= max_order; n++) {
            const Poly Hn = K_current.homogeneous(n);
            if (Hn.empty()) continue;

            Poly Gn = compute_generator(Hn, nu, denom_tol);
            res.G[n] = Gn;

            K_current = apply_lie_transform_H(K_current, Gn, max_order, n);
            res.K_normal.add_poly(extract_resonant_terms(K_current, nu, n, resonance_tol));
        }

        return res;
    }

    ActionAngle evaluate_action_angle(const std::array<Cplx, 6>& qp) {
        const Cplx& q1 = qp[0];
        const Cplx& q2 = qp[1];
        const Cplx& q3 = qp[2];
        const Cplx& p1 = qp[3];
        const Cplx& p2 = qp[4];
        const Cplx& p3 = qp[5];

        ActionAngle out{};
        out.I[0] = q1 * p1;
        out.I[1] = q2 * (Cplx{ 0, 1 } *p2);  // q2 * i p2
        out.I[2] = q3 * (Cplx{ 0, 1 } *p3);

        // phi1 = Re(0.5 * log(q1/p1))
        const Cplx ratio1 = q1 / p1;
        const Cplx log_ratio = std::log(ratio1);
        out.phi[0] = Real(0.5L) * log_ratio.real();

        out.phi[1] = std::arg(q2);
        out.phi[2] = std::arg(q3);
        return out;
    }


    ActionsExpansionFull extract_actions_expansion_full(const Poly& K, int max_I_degree) {
        ActionsExpansionFull out;

        auto add_Iterm = [&](int a1, int a2, int a3, Cplx c) {
            if (c == Cplx{ 0,0 }) return;
            int d = a1 + a2 + a3;
            if (d == 0 || d > max_I_degree) return;
            IExp3 key{ a1,a2,a3 };
            auto it = out.terms.find(key);
            if (it == out.terms.end()) out.terms.emplace(key, c);
            else {
                it->second += c;
                if (it->second == Cplx{ 0,0 }) out.terms.erase(it);
            }
            };

        for (const auto& kv : K.terms()) {
            const Exp6& e = kv.first;
            Cplx c = kv.second;

            int e1 = e[0], e2 = e[1], e3 = e[2];
            int e4 = e[3], e5 = e[4], e6 = e[5];

            // 只能处理 q1^n p1^n / q2^n p2^n / q3^n p3^n 的“配对”情况，
            // 如果指数不匹配，说明这项不是纯 I 的函数（在正规形里这类项应当已被消掉）
            if (e1 != e4 || e2 != e5 || e3 != e6) continue;

            int a1 = e1;   // I1^a1
            int a2 = e2;   // I2^a2
            int a3 = e3;   // I3^a3

            // I2 = i q2 p2, I3 = i q3 p3
            // => (q2 p2)^{a2}(q3 p3)^{a3} = (-i)^{a2+a3} I2^{a2} I3^{a3}
            int k = a2 + a3; // 幂次
            Cplx factor;
            switch ((k % 4 + 4) % 4) {
            case 0: factor = Cplx{ 1,0 }; break;
            case 1: factor = Cplx{ 0,-1 }; break;   // -i
            case 2: factor = Cplx{ -1,0 }; break;
            case 3: factor = Cplx{ 0,1 }; break;    // +i
            }
            Cplx cI = c * factor;

            add_Iterm(a1, a2, a3, cI);
        }

        return out;
    }


    ActionsExpansion extract_actions_expansion(const Poly& K) {
        ActionsExpansion a{};

        auto get = [&](std::uint8_t e0, std::uint8_t e1, std::uint8_t e2,
            std::uint8_t e3, std::uint8_t e4, std::uint8_t e5) -> Cplx {
                Exp6 key{ e0,e1,e2,e3,e4,e5 };
                auto it = K.terms().find(key);
                if (it == K.terms().end()) return Cplx{ 0,0 };
                return it->second;
            };

        // 线性项
        Cplx c_q1p1 = get(1, 0, 0, 1, 0, 0);
        Cplx c_q2p2 = get(0, 1, 0, 0, 1, 0);
        Cplx c_q3p3 = get(0, 0, 1, 0, 0, 1);

        a.cI1 = c_q1p1;
        a.cI2 = -Cplx{ 0,1 } *c_q2p2; // 因为 q2 p2 = -i I2 => c_I2 = -i*c_q2p2
        a.cI3 = -Cplx{ 0,1 } *c_q3p3;

        // 二次自项
        Cplx c_q1q1p1p1 = get(2, 0, 0, 2, 0, 0);
        Cplx c_q2q2p2p2 = get(0, 2, 0, 0, 2, 0);
        Cplx c_q3q3p3p3 = get(0, 0, 2, 0, 0, 2);

        a.cI1I1 = c_q1q1p1p1;
        a.cI2I2 = -c_q2q2p2p2;  // I2^2 = -q2^2 p2^2
        a.cI3I3 = -c_q3q3p3p3;

        // 交叉项
        Cplx c_q1q2p1p2 = get(1, 1, 0, 1, 1, 0);
        Cplx c_q1q3p1p3 = get(1, 0, 1, 1, 0, 1);
        Cplx c_q2q3p2p3 = get(0, 1, 1, 0, 1, 1);

        // I1 I2 = i q1q2p1p2 => q1q2p1p2 = -i I1I2
        a.cI1I2 = -Cplx{ 0,1 } *c_q1q2p1p2;
        // I1 I3 = i q1q3p1p3
        a.cI1I3 = -Cplx{ 0,1 } *c_q1q3p1p3;
        // I2 I3 = - q2q3p2p3
        a.cI2I3 = -c_q2q3p2p3;

        return a;
    }

    // 对单个坐标 poly 做 exp(L_G) 变换，和 apply_lie_transform_H 一样的阶数策略
    Poly apply_lie_transform_coord(const Poly& coord,
        const Poly& G, int max_order, int n) {
        Poly out = coord;
        Poly Lprev = poisson_bracket(coord, G).truncated(max_order); // L_G coord
        int k = 1;
        Real factorial = 1.0L;

        while (!Lprev.empty()) {
            int min_deg = 1 + k * (n - 2);        // 坐标最低是一阶
            if (min_deg > max_order) break;

            factorial *= static_cast<Real>(k);
            out.add_poly(Lprev, Cplx{ 1.0L / factorial, 0 });

            Poly Lnext = poisson_bracket(Lprev, G).truncated(max_order);
            Lprev = std::move(Lnext);
            ++k;
        }

        return out.truncated(max_order);
    }

    // forward: (q,p)_complex -> (Q,P)_normal_form
    void build_coordinate_transforms_forward(
        const nf::NormalFormResult& res, int max_order,
        nf::Poly& Q1, nf::Poly& Q2, nf::Poly& Q3,
        nf::Poly& P1, nf::Poly& P2, nf::Poly& P3) {

        using namespace nf;
        // 1. 初始化为恒等变换
        Q1 = Poly{}; Q1.add_term(Exp6{ 1,0,0,0,0,0 }, Cplx{ 1,0 }); // q1
        Q2 = Poly{}; Q2.add_term(Exp6{ 0,1,0,0,0,0 }, Cplx{ 1,0 }); // q2
        Q3 = Poly{}; Q3.add_term(Exp6{ 0,0,1,0,0,0 }, Cplx{ 1,0 }); // q3
        P1 = Poly{}; P1.add_term(Exp6{ 0,0,0,1,0,0 }, Cplx{ 1,0 }); // p1
        P2 = Poly{}; P2.add_term(Exp6{ 0,0,0,0,1,0 }, Cplx{ 1,0 }); // p2
        P3 = Poly{}; P3.add_term(Exp6{ 0,0,0,0,0,1 }, Cplx{ 1,0 }); // p3

        // 2. 按阶数依次应用每个 G_n
        for (int n = 3; n < (int)res.G.size(); ++n) {
            const Poly& Gn = res.G[n];
            if (Gn.empty()) continue;

            Q1 = apply_lie_transform_coord(Q1, Gn, max_order, n);
            Q2 = apply_lie_transform_coord(Q2, Gn, max_order, n);
            Q3 = apply_lie_transform_coord(Q3, Gn, max_order, n);
            P1 = apply_lie_transform_coord(P1, Gn, max_order, n);
            P2 = apply_lie_transform_coord(P2, Gn, max_order, n);
            P3 = apply_lie_transform_coord(P3, Gn, max_order, n);
        }
    }

}  // namespace nf

