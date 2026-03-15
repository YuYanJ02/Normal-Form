#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>

namespace nf {

    using Real = long double;
    using Cplx = std::complex<Real>;

    // Exponent key for monomial: q1,q2,q3,p1,p2,p3
    using Exp6 = std::array<std::uint8_t, 6>;

    struct Exp6Hash {
        std::size_t operator()(const Exp6& e) const noexcept;
    };

    inline int total_degree(const Exp6& e) {
        return int(e[0]) + int(e[1]) + int(e[2]) + int(e[3]) + int(e[4]) + int(e[5]);
    }

    class Poly {
    public:
        using Map = std::unordered_map<Exp6, Cplx, Exp6Hash>;

        Poly() = default;
        explicit Poly(Cplx c);  // constant

        const Map& terms() const { return terms_; }
        Map& terms() { return terms_; }

        bool empty() const { return terms_.empty(); }

        void add_term(const Exp6& e, Cplx c);
        void add_poly(const Poly& other, Cplx scale = Cplx{ 1, 0 });

        Poly scaled(Cplx s) const;
        Poly truncated(int max_degree) const;
        Poly homogeneous(int degree) const;

        void prune(Real abs_threshold);

        std::string to_string(const std::array<std::string, 6>& var_names,
            int max_terms = 50) const;

    private:
        Map terms_;
    };

    // Canonical Poisson bracket {f,g} on (q1,q2,q3,p1,p2,p3).
    Poly poisson_bracket(const Poly& f, const Poly& g);

    struct NormalFormResult {
        Poly K_normal;
        std::vector<Poly> G;  // G[n] valid for n=0..N (G[0],G[1],G[2] empty)
    };

    // Compute G_n from homogeneous part H_n and nu = [lambda, i*omega_y, i*omega_z].
    Poly compute_generator(const Poly& Hn, const std::array<Cplx, 3>& nu,
        Real denom_tol = 1e-18L);

    // Extract resonant terms in homogeneous order n, using tolerance on |<kp-kq,nu>|.
    Poly extract_resonant_terms(const Poly& K_current, const std::array<Cplx, 3>& nu,
        int order_n, Real resonance_tol = 1e-18L);

    // Apply exp(L_G) to K truncated to max_order (up to cubic in Lie series).
    Poly apply_lie_transform_H(const Poly& K, const Poly& G, int max_order, int n);

    NormalFormResult lie_transformation(const Poly& H_complex,
        const std::array<Cplx, 3>& nu,
        int max_order,
        Real denom_tol = 1e-18L,
        Real resonance_tol = 1e-18L);

    // Helpers for action-angle evaluation from normal-form complex coords.
    // qp = [q1,q2,q3,p1,p2,p3] (q2,q3,p2,p3 complex generally)
    struct ActionAngle {
        std::array<Cplx, 3> I;
        std::array<Real, 3> phi;  // phi1 real-valued (we return Re(0.5*log(q1/p1)))
    };

    ActionAngle evaluate_action_angle(const std::array<Cplx, 6>& qp);

    Poly load_poly_from_file(const std::string& path);

    struct ActionsExpansion {
        Cplx cI1, cI2, cI3;
        Cplx cI1I1, cI2I2, cI3I3;
        Cplx cI1I2, cI1I3, cI2I3;
    };

    ActionsExpansion extract_actions_expansion(const Poly& K);

    using IExp3 = std::array<int, 3>;

    struct IExp3Hash {
        std::size_t operator()(const IExp3& e) const noexcept {
            std::size_t h = 1469598103934665603ull;
            for (int v : e) {
                h ^= static_cast<std::size_t>(v);
                h *= 1099511628211ull;
            }
            return h;
        }
    };

    struct ActionsExpansionFull {
        std::unordered_map<IExp3, Cplx, IExp3Hash> terms;
    };

    ActionsExpansionFull extract_actions_expansion_full(const Poly& K, int max_I_degree);

    // 坐标 Lie 变换：复数坐标 (q,p) -> 正规形坐标 (Q,P)
    Poly apply_lie_transform_coord(const Poly& coord,
        const Poly& G,
        int max_order,
        int n);

    void build_coordinate_transforms_forward(
        const NormalFormResult& res,
        int max_order,
        Poly& Q1, Poly& Q2, Poly& Q3,
        Poly& P1, Poly& P2, Poly& P3);

}  // namespace nf

