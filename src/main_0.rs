// first тензор метрики и обращение
 
// use std::f64::consts::PI;
// use autodiff::Float;
// use nalgebra::{Matrix4, Vector4};
// use ndarray::Array2;

// fn main() {
//     // Параметры модели
//     let r = 1.0;
//     let theta = PI / 4.0;
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0; // натуральные единицы

//     // Компоненты метрики
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6) / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_rtheta = -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);
//     let g_zz = 1.0;

//     // Метрика (4x4) в координатах (t, r, θ, z)
//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(1, 2)] = g_rtheta;
//     g[(2, 1)] = g_rtheta;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = g_zz;

//     // Вывод метрики
//     println!("Metric tensor g:");
//     println!("{:.5}", g);

//     // Обратная метрика
//     let g_inv = g.try_inverse().expect("Metric is not invertible");

//     println!("\nInverse metric g^(-1):");
//     println!("{:.5}", g_inv);

//     // TODO: здесь можно будет добавить вычисление символов Кристоффеля,
//     // тензора Риччи и тензора Эйнштейна G_{μν}
// }


//second для проверки тензора Эйнштейна

// use std::f64::consts::PI;
// use nalgebra::{Matrix4, Vector4};
// use ndarray::Array3;

// fn g_components(r: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6) / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_rtheta = -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);
//     let g_zz = 1.0;

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(1, 2)] = g_rtheta;
//     g[(2, 1)] = g_rtheta;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = g_zz;
//     g
// }

// // Производные метрики по r по компонентам: ∂_r g_{μν}
// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> [Matrix4<f64>; 1] {
//     let g_plus = g_components(r + h, r0, n, c);
//     let g_minus = g_components(r - h, r0, n, c);
//     let mut dg_dr = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg_dr[(mu, nu)] = (g_plus[(mu, nu)] - g_minus[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     [dg_dr]
// }

// fn christoffel_symbols(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &[Matrix4<f64>]) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     // ∂_μ g_{νσ} = dg[0] (только ∂_r, index 1)
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let d1 = if mu == 1 { dg[0][(nu, sigma)] } else { 0.0 };
//                     let d2 = if nu == 1 { dg[0][(mu, sigma)] } else { 0.0 };
//                     let d3 = if sigma == 1 { dg[0][(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (d1 + d2 - d3);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// fn ricci_tensor(gamma: &Array3<f64>, h: f64) -> Matrix4<f64> {
//     let mut ricci = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             let mut sum = 0.0;
//             for lambda in 0..4 {
//                 for rho in 0..4 {
//                     // Только производные по r (index 1)
//                     let d_gamma_l_mn = (gamma[[lambda, mu, nu]] - gamma[[lambda, mu, nu]]) / h; // = 0
//                     let d_gamma_l_ml = (gamma[[lambda, mu, lambda]] - gamma[[lambda, mu, lambda]]) / h; // = 0
//                     let term1 = d_gamma_l_mn - d_gamma_l_ml;
//                     let term2 = gamma[[lambda, lambda, rho]] * gamma[[rho, mu, nu]];
//                     let term3 = gamma[[lambda, nu, rho]] * gamma[[rho, mu, lambda]];
//                     sum += term1 + term2 - term3;
//                 }
//             }
//             ricci[(mu, nu)] = sum;
//         }
//     }
//     ricci
// }

// fn main() {
//     let r = 1.0;
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     let g = g_components(r, r0, n, c);
//     let g_inv = g.try_inverse().unwrap();
//     let dg = metric_derivatives(r, r0, n, c, h);
//     let gamma = christoffel_symbols(&g, &g_inv, &dg);
//     let ricci = ricci_tensor(&gamma, h);

//     // Скалярная кривизна R
//     let mut scalar_r = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             scalar_r += g_inv[(mu, nu)] * ricci[(mu, nu)];
//         }
//     }

//     // Тензор Эйнштейна G = R_{μν} - ½ R g_{μν}
//     let mut einstein = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             einstein[(mu, nu)] = ricci[(mu, nu)] - 0.5 * scalar_r * g[(mu, nu)];
//         }
//     }

//     println!("Einstein tensor G_μν:");
//     println!("{:.5}", einstein);
// }



// #3 отклонения после проверки тензора Эйнштейна
// use std::f64::consts::PI;
// use nalgebra::{Matrix4};
// use ndarray::Array3;
// use plotters::prelude::*; // для графика

// // Метрика без g_{rθ}, как раньше
// fn g_components(r: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6) / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);
//     let g_zz = 1.0;

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = g_zz;
//     g
// }

// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> [Matrix4<f64>; 1] {
//     let g_plus = g_components(r + h, r0, n, c);
//     let g_minus = g_components(r - h, r0, n, c);
//     let mut dg_dr = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg_dr[(mu, nu)] = (g_plus[(mu, nu)] - g_minus[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     [dg_dr]
// }

// fn christoffel_symbols(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &[Matrix4<f64>]) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let d1 = if mu == 1 { dg[0][(nu, sigma)] } else { 0.0 };
//                     let d2 = if nu == 1 { dg[0][(mu, sigma)] } else { 0.0 };
//                     let d3 = if sigma == 1 { dg[0][(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (d1 + d2 - d3);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// fn ricci_tensor(gamma: &Array3<f64>, h: f64) -> Matrix4<f64> {
//     let mut ricci = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             let mut sum = 0.0;
//             for lambda in 0..4 {
//                 for rho in 0..4 {
//                     let term1 = 0.0; // dΓ - dΓ = 0 (по r только)
//                     let term2 = gamma[[lambda, lambda, rho]] * gamma[[rho, mu, nu]];
//                     let term3 = gamma[[lambda, nu, rho]] * gamma[[rho, mu, lambda]];
//                     sum += term1 + term2 - term3;
//                 }
//             }
//             ricci[(mu, nu)] = sum;
//         }
//     }
//     ricci
// }

// fn g_tt(r: f64) -> f64 {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     let g = g_components(r, r0, n, c);
//     let g_inv = g.try_inverse().unwrap();
//     let dg = metric_derivatives(r, r0, n, c, h);
//     let gamma = christoffel_symbols(&g, &g_inv, &dg);
//     let ricci = ricci_tensor(&gamma, h);

//     let mut scalar_r = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             scalar_r += g_inv[(mu, nu)] * ricci[(mu, nu)];
//         }
//     }

//     let g_tt = ricci[(0, 0)] - 0.5 * scalar_r * g[(0, 0)];
//     g_tt
// }

// fn main() -> Result<(), Box<dyn std::error::Error>> {
//     let root = BitMapBackend::new("gtt_plot.png", (800, 600)).into_drawing_area();
//     root.fill(&WHITE)?;
//     let mut chart = ChartBuilder::on(&root)
//         .caption("G_tt(r)", ("sans-serif", 30))
//         .margin(10)
//         .x_label_area_size(40)
//         .y_label_area_size(50)
//         .build_cartesian_2d(0.5..10.0, -1.0..0.2)?;

//     chart.configure_mesh().draw()?;

//     chart.draw_series(LineSeries::new(
//         (50..200).map(|i| {
//             let r = i as f64 / 15.0;
//             (r, g_tt(r))
//         }),
//         &RED,
//     ))?;

//     Ok(())
// }


// #4 геодезические уравнения HyperTwist

// use std::f64::consts::PI;
// use nalgebra::{Vector4, Matrix4};
// use ndarray::Array3;

// /// Типы данных
// type State = Vector4<f64>;   // координаты (t, r, θ, z)
// type DState = Vector4<f64>;  // 4-скорость

// /// Метрика HyperTwist
// fn metric(r: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_rtheta = -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(1, 2)] = g_rtheta;
//     g[(2, 1)] = g_rtheta;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// /// Производные ∂_r g_{μν}
// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> [Matrix4<f64>; 1] {
//     let gp = metric(r + h, r0, n, c);
//     let gm = metric(r - h, r0, n, c);
//     let mut dg = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg[(mu, nu)] = (gp[(mu, nu)] - gm[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     [dg]
// }

// /// Символы Кристоффеля Γ^λ_{μν}
// fn christoffel(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &[Matrix4<f64>]) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let dmu = if mu == 1 { dg[0][(nu, sigma)] } else { 0.0 };
//                     let dnu = if nu == 1 { dg[0][(mu, sigma)] } else { 0.0 };
//                     let dsigma = if sigma == 1 { dg[0][(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// /// Геодезическое ускорение d²x^μ/dτ² = -Γ^μ_{νρ} dx^ν dx^ρ
// fn geodesic_rhs(dx: &DState, gamma: &Array3<f64>) -> DState {
//     let mut ddx = Vector4::<f64>::zeros();
//     for mu in 0..4 {
//         let mut acc = 0.0;
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 acc -= gamma[[mu, nu, rho]] * dx[nu] * dx[rho];
//             }
//         }
//         ddx[mu] = acc;
//     }
//     ddx
// }

// /// main()
// fn main() {
//     let r = 1.0;
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     let g = metric(r, r0, n, c);
//     let g_inv = g.try_inverse().unwrap();
//     let dg = metric_derivatives(r, r0, n, c, h);
//     let gamma = christoffel(&g, &g_inv, &dg);

//     // Пример 4-скорости: (dt/dτ, dr/dτ, dθ/dτ, dz/dτ)
//     let dx = Vector4::new(1.0, 0.0, 0.5, 0.0);

//     let ddx = geodesic_rhs(&dx, &gamma);

//     println!("Геодезическое ускорение d²x^μ/dτ²:");
//     for mu in 0..4 {
//         println!("  μ = {} → {}", mu, ddx[mu]);
//     }
// }


// #5  ✅ Rust-код: численно посчитать v(r) из геодезии
// use std::f64::consts::PI;
// use nalgebra::{Matrix4};

// fn metric_components(r: f64, r0: f64, n: f64, c: f64) -> (f64, f64, f64, f64) {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_rtheta = -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);
//     (f, g_rr, g_rtheta, g_thetatheta)
// }

// fn gamma_r_tt(r: f64, r0: f64, n: f64, c: f64, h: f64) -> f64 {
//     let (f1, _, _, _) = metric_components(r + h, r0, n, c);
//     let (f2, _, _, _) = metric_components(r - h, r0, n, c);
//     let df = (f1 - f2) / (2.0 * h);
//     let g_rr = metric_components(r, r0, n, c).1;
//     0.5 * df * c * c / g_rr
// }

// fn gamma_r_thetatheta(r: f64, r0: f64, n: f64, c: f64, h: f64) -> f64 {
//     let (_, _, _, g_th1) = metric_components(r + h, r0, n, c);
//     let (_, _, _, g_th2) = metric_components(r - h, r0, n, c);
//     let dg = (g_th1 - g_th2) / (2.0 * h);
//     let g_rr = metric_components(r, r0, n, c).1;
//     -0.5 * dg / g_rr
// }

// fn orbital_velocity(r: f64, r0: f64, n: f64, c: f64, h: f64) -> f64 {
//     let grtt = gamma_r_tt(r, r0, n, c, h);
//     let grthth = gamma_r_thetatheta(r, r0, n, c, h);
//     let ratio = grthth / grtt;
//     // println!("r = {:.2} | Γ^r_tt = {:.5}, Γ^r_thth = {:.5}, ratio = {:.5}", r, grtt, grthth, ratio);

//     if ratio > 0.0 {
//         (r * (ratio).sqrt())
//     } else {
//         0.0 // вне допустимой области
//     }
// }

// fn main() {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     println!("r\tv(r)");
//     for i in 1..100 {
//         let r = i as f64 * 0.1;
//         let v = orbital_velocity(r, r0, n, c, h);
//         println!("{:.2}\t{:.5}", r, v);
//     }
// }


// #6 - тензор ричи

// use std::f64::consts::PI;
// use nalgebra::{Matrix4};
// use ndarray::Array3;

// /// Метрика HyperTwist (упрощённая, без g_{rθ})
// fn metric(r: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// /// Производные по r
// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> Matrix4<f64> {
//     let gp = metric(r + h, r0, n, c);
//     let gm = metric(r - h, r0, n, c);
//     let mut dg = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg[(mu, nu)] = (gp[(mu, nu)] - gm[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     dg
// }

// /// Символы Кристоффеля: только ∂_r g_{μν}
// fn christoffel(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &Matrix4<f64>) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let dmu = if mu == 1 { dg[(nu, sigma)] } else { 0.0 };
//                     let dnu = if nu == 1 { dg[(mu, sigma)] } else { 0.0 };
//                     let dsigma = if sigma == 1 { dg[(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// /// Тензор Риччи R_{μν}
// fn ricci_tensor(gamma: &Array3<f64>, h: f64) -> Matrix4<f64> {
//     let mut ricci = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             let mut sum = 0.0;
//             for lambda in 0..4 {
//                 for rho in 0..4 {
//                     let term1 = 0.0; // производные не считаем, но можно добавить
//                     let term2 = gamma[[lambda, lambda, rho]] * gamma[[rho, mu, nu]];
//                     let term3 = gamma[[lambda, nu, rho]] * gamma[[rho, mu, lambda]];
//                     sum += term1 + term2 - term3;
//                 }
//             }
//             ricci[(mu, nu)] = sum;
//         }
//     }
//     ricci
// }

// /// Скалярная кривизна R = g^{μν} R_{μν}
// fn scalar_curvature(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
//     let mut r = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             r += g_inv[(mu, nu)] * ricci[(mu, nu)];
//         }
//     }
//     r
// }

// /// Инвариант Риччи R_{μν} R^{μν}
// fn ricci_invariant(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
//     let mut sum = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             for alpha in 0..4 {
//                 for beta in 0..4 {
//                     sum += g_inv[(mu, alpha)] * g_inv[(nu, beta)] * ricci[(mu, nu)] * ricci[(alpha, beta)];
//                 }
//             }
//         }
//     }
//     sum
// }

// fn main() {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     println!("r\tR(r)\t\tRicci^2");

//     for i in 1..100 {
//         let r = i as f64 * 0.1;

//         let g = metric(r, r0, n, c);
//         let g_inv = g.try_inverse().unwrap();
//         let dg = metric_derivatives(r, r0, n, c, h);
//         let gamma = christoffel(&g, &g_inv, &dg);
//         let ricci = ricci_tensor(&gamma, h);

//         let scalar_r = scalar_curvature(&g_inv, &ricci);
//         let ricci_sq = ricci_invariant(&g_inv, &ricci);

//         println!("{:.2}\t{:.6}\t{:.6}", r, scalar_r, ricci_sq);
//     }
// }

// #7 - гравитационное линзирование

// use std::f64::consts::PI;
// use std::fs::File;
// use std::io::Write;
// use nalgebra::{Vector4, Matrix4};
// use ndarray::Array3;

// type State = Vector4<f64>; // (t, r, θ, z)
// type DState = Vector4<f64>; // d/dλ (t, r, θ, z)

// fn metric(r: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_rtheta = -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(1, 2)] = g_rtheta;
//     g[(2, 1)] = g_rtheta;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> [Matrix4<f64>; 1] {
//     let gp = metric(r + h, r0, n, c);
//     let gm = metric(r - h, r0, n, c);
//     let mut dg = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg[(mu, nu)] = (gp[(mu, nu)] - gm[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     [dg]
// }

// fn christoffel(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &[Matrix4<f64>]) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let dmu = if mu == 1 { dg[0][(nu, sigma)] } else { 0.0 };
//                     let dnu = if nu == 1 { dg[0][(mu, sigma)] } else { 0.0 };
//                     let dsigma = if sigma == 1 { dg[0][(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// fn geodesic_rhs(dx: &DState, gamma: &Array3<f64>) -> DState {
//     let mut ddx = Vector4::<f64>::zeros();
//     for mu in 0..4 {
//         let mut acc = 0.0;
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 acc -= gamma[[mu, nu, rho]] * dx[nu] * dx[rho];
//             }
//         }
//         ddx[mu] = acc;
//     }
//     ddx
// }

// fn rk4_step(
//     x: &mut State,
//     dx: &mut DState,
//     r0: f64,
//     n: f64,
//     c: f64,
//     h: f64,
//     dλ: f64,
// ) {
//     let r = x[1].max(1e-5); // избегаем r=0
//     let g = metric(r, r0, n, c);
//     let g_inv = g.try_inverse().unwrap();
//     let dg = metric_derivatives(r, r0, n, c, h);
//     let gamma = christoffel(&g, &g_inv, &dg);

//     let k1 = geodesic_rhs(dx, &gamma);
//     let dx2 = *dx + 0.5 * dλ * k1;
//     let k2 = geodesic_rhs(&dx2, &gamma);
//     let dx3 = *dx + 0.5 * dλ * k2;
//     let k3 = geodesic_rhs(&dx3, &gamma);
//     let dx4 = *dx + dλ * k3;
//     let k4 = geodesic_rhs(&dx4, &gamma);

//     *dx += dλ * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
//     *x += dλ * *dx;
// }

// fn main() -> std::io::Result<()> {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;
//     let dλ = 0.01;

//     // Фотон стартует слева от центра, летит вправо
//     let mut x = Vector4::new(0.0, 5.0, PI, 0.0);       // t, r, θ, z
//     let mut dx = Vector4::new(1.0, -1.0, 0.3, 0.0);     // светоподобная 4-скорость

//     // Нормировка: g_{μν} dx^μ dx^ν = 0
//     {
//         let g = metric(x[1], r0, n, c);
//         let gdot = dx.transpose() * g * dx;
//         let norm = gdot[(0, 0)];
//         dx /= norm.sqrt(); // масштабируем до светоподобной траектории
//     }

//     let mut file = File::create("lens_trajectory.csv")?;
//     writeln!(file, "x,y")?;

//     for _ in 0..3000 {
//         let r = x[1];
//         let theta = x[2];
//         let px = r * theta.cos();
//         let py = r * theta.sin();
//         writeln!(file, "{},{}", px, py)?;
//         rk4_step(&mut x, &mut dx, r0, n, c, h, dλ);

//         // остановка после пролёта
//         if px > 5.0 && r > 5.0 {
//             break;
//         }
//     }

//     println!("✅ Траектория фотона записана в lens_trajectory.csv");
//     Ok(())
// }

// #8 - анализ угла отклонения света
// use std::fs::File;
// use std::io::{BufReader, BufRead};
// use std::f64::consts::PI;

// fn angular_difference(theta1: f64, theta2: f64) -> f64 {
//     let mut d = theta2 - theta1;
//     while d > PI { d -= 2.0 * PI; }
//     while d < -PI { d += 2.0 * PI; }
//     d
// }

// fn main() -> std::io::Result<()> {
//     let file = File::open("lens_trajectory.csv")?;
//     let reader = BufReader::new(file);
//     let lines: Vec<_> = reader.lines().skip(1)
//         .filter_map(Result::ok)
//         .collect();

//     if lines.len() < 2 {
//         println!("Недостаточно точек для анализа.");
//         return Ok(());
//     }

//     let (x1, y1) = {
//         let parts: Vec<_> = lines[0].split(',').collect();
//         (parts[0].parse::<f64>().unwrap(), parts[1].parse::<f64>().unwrap())
//     };

//     // Найдём последнюю НЕ-NaN строку
//     let (x2, y2) = lines.iter().rev()
//         .filter_map(|line| {
//             let parts: Vec<_> = line.split(',').collect();
//             let x = parts[0].parse::<f64>().ok()?;
//             let y = parts[1].parse::<f64>().ok()?;
//             if x.is_finite() && y.is_finite() {
//                 Some((x, y))
//             } else {
//                 None
//             }
//         })
//         .next()
//         .expect("Нет допустимых финальных точек");

//     let theta_in = y1.atan2(x1);
//     let theta_out = y2.atan2(x2);
//     let delta = angular_difference(theta_in, theta_out).abs().to_degrees();

//     println!("🔭 Угол входа  θ_in  = {:.4} рад", theta_in);
//     println!("🔭 Угол выхода θ_out = {:.4} рад", theta_out);
//     println!("➡️  Угол отклонения света Δθ = {:.6}°", delta);

//     Ok(())
// }

// #9 - тензор римана 

// use std::f64::consts::PI;
// use nalgebra::{Matrix4};
// use ndarray::{Array3, Array4};

// fn metric(r: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> Matrix4<f64> {
//     let gp = metric(r + h, r0, n, c);
//     let gm = metric(r - h, r0, n, c);
//     let mut dg = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg[(mu, nu)] = (gp[(mu, nu)] - gm[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     dg
// }

// fn christoffel(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &Matrix4<f64>) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let dmu = if mu == 1 { dg[(nu, sigma)] } else { 0.0 };
//                     let dnu = if nu == 1 { dg[(mu, sigma)] } else { 0.0 };
//                     let dsigma = if sigma == 1 { dg[(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// fn ricci_tensor(gamma: &Array3<f64>, h: f64) -> Matrix4<f64> {
//     let mut ricci = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             let mut sum = 0.0;
//             for lambda in 0..4 {
//                 for rho in 0..4 {
//                     let term1 = 0.0;
//                     let term2 = gamma[[lambda, lambda, rho]] * gamma[[rho, mu, nu]];
//                     let term3 = gamma[[lambda, nu, rho]] * gamma[[rho, mu, lambda]];
//                     sum += term1 + term2 - term3;
//                 }
//             }
//             ricci[(mu, nu)] = sum;
//         }
//     }
//     ricci
// }

// fn scalar_curvature(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
//     let mut r = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             r += g_inv[(mu, nu)] * ricci[(mu, nu)];
//         }
//     }
//     r
// }

// fn ricci_invariant(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
//     let mut sum = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             for alpha in 0..4 {
//                 for beta in 0..4 {
//                     sum += g_inv[(mu, alpha)] * g_inv[(nu, beta)] * ricci[(mu, nu)] * ricci[(alpha, beta)];
//                 }
//             }
//         }
//     }
//     sum
// }

// fn riemann_tensor(gamma: &Array3<f64>, h: f64) -> Array4<f64> {
//     let mut riemann = Array4::<f64>::zeros((4, 4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     let d_gamma_rho = 0.0; // численные производные можно уточнить
//                     let d_gamma_sigma = 0.0;
//                     let mut term3 = 0.0;
//                     let mut term4 = 0.0;
//                     for nu in 0..4 {
//                         term3 += gamma[[lambda, nu, rho]] * gamma[[nu, mu, sigma]];
//                         term4 += gamma[[lambda, nu, sigma]] * gamma[[nu, mu, rho]];
//                     }
//                     riemann[[lambda, mu, rho, sigma]] =
//                         d_gamma_rho - d_gamma_sigma + term3 - term4;
//                 }
//             }
//         }
//     }
//     riemann
// }

// fn kretschmann_invariant(g_inv: &Matrix4<f64>, riemann: &Array4<f64>) -> f64 {
//     let mut sum = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     for alpha in 0..4 {
//                         for beta in 0..4 {
//                             for gamma in 0..4 {
//                                 for delta in 0..4 {
//                                     sum += g_inv[(mu, alpha)] * g_inv[(nu, beta)] *
//                                             g_inv[(rho, gamma)] * g_inv[(sigma, delta)] *
//                                             riemann[[mu, nu, rho, sigma]] * riemann[[alpha, beta, gamma, delta]];
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     sum
// }

// fn main() {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     println!("r\tR\t\tRicci^2\t\tKretsch");

//     for i in 1..100 {
//         let r = i as f64 * 0.1;

//         let g = metric(r, r0, n, c);
//         let g_inv = g.try_inverse().unwrap();
//         let dg = metric_derivatives(r, r0, n, c, h);
//         let gamma = christoffel(&g, &g_inv, &dg);
//         let ricci = ricci_tensor(&gamma, h);
//         let scalar_r = scalar_curvature(&g_inv, &ricci);
//         let ricci_sq = ricci_invariant(&g_inv, &ricci);
//         let riemann = riemann_tensor(&gamma, h);
//         let kretsch = kretschmann_invariant(&g_inv, &riemann);

//         println!("{:.2}\t{:.6}\t{:.6}\t{:.6}", r, scalar_r, ricci_sq, kretsch);
//     }
// }


// # 10 - метрика действительно не зависит от угла
// use std::f64::consts::PI;
// use nalgebra::{Matrix4};
// use ndarray::{Array3, Array4};

// fn metric(r: f64, _theta: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// fn is_axially_symmetric(r: f64, r0: f64, n: f64, c: f64, epsilon: f64) -> bool {
//     let theta1 = 0.0;
//     let theta2 = PI / 2.0;
//     let g1 = metric(r, theta1, r0, n, c);
//     let g2 = metric(r, theta2, r0, n, c);

//     for mu in 0..4 {
//         for nu in 0..4 {
//             if (g1[(mu, nu)] - g2[(mu, nu)]).abs() > epsilon {
//                 println!("❌ Asymmetry at r = {:.2}: g[{},{}] differs by {:.2e}", r, mu, nu, (g1[(mu, nu)] - g2[(mu, nu)]));
//                 return false;
//             }
//         }
//     }
//     true
// }

// fn main() {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let epsilon = 1e-10;

//     println!("\n=== Axial symmetry check ===");
//     for i in 1..100 {
//         let r = i as f64 * 0.1;
//         if !is_axially_symmetric(r, r0, n, c, epsilon) {
//             println!("❌ Symmetry broken at r = {:.2}", r);
//             return;
//         }
//     }
//     println!("✅ Axial symmetry confirmed for all r in [0.1, 10.0]");
// }


// #11 - тензор Вейля

// use std::f64::consts::PI;
// use nalgebra::Matrix4;
// use ndarray::{Array3, Array4};

// fn metric(r: f64, _theta: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> Matrix4<f64> {
//     let gp = metric(r + h, 0.0, r0, n, c);
//     let gm = metric(r - h, 0.0, r0, n, c);
//     let mut dg = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg[(mu, nu)] = (gp[(mu, nu)] - gm[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     dg
// }

// fn christoffel(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &Matrix4<f64>) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let dmu = if mu == 1 { dg[(nu, sigma)] } else { 0.0 };
//                     let dnu = if nu == 1 { dg[(mu, sigma)] } else { 0.0 };
//                     let dsigma = if sigma == 1 { dg[(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// fn riemann_tensor(gamma: &Array3<f64>) -> Array4<f64> {
//     let mut riemann = Array4::<f64>::zeros((4, 4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     let mut term3 = 0.0;
//                     let mut term4 = 0.0;
//                     for nu in 0..4 {
//                         term3 += gamma[[lambda, nu, rho]] * gamma[[nu, mu, sigma]];
//                         term4 += gamma[[lambda, nu, sigma]] * gamma[[nu, mu, rho]];
//                     }
//                     riemann[[lambda, mu, rho, sigma]] = term3 - term4;
//                 }
//             }
//         }
//     }
//     riemann
// }

// fn ricci_tensor(gamma: &Array3<f64>) -> Matrix4<f64> {
//     let mut ricci = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             let mut sum = 0.0;
//             for lambda in 0..4 {
//                 sum += gamma[[lambda, mu, nu]] * gamma[[lambda, nu, mu]]
//                     - gamma[[lambda, lambda, mu]] * gamma[[nu, nu, lambda]];
//             }
//             ricci[(mu, nu)] = sum;
//         }
//     }
//     ricci
// }

// fn scalar_curvature(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
//     let mut r = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             r += g_inv[(mu, nu)] * ricci[(mu, nu)];
//         }
//     }
//     r
// }

// fn weyl_tensor(
//     g: &Matrix4<f64>,
//     g_inv: &Matrix4<f64>,
//     ricci: &Matrix4<f64>,
//     scalar_r: f64,
//     riemann: &Array4<f64>,
// ) -> Array4<f64> {
//     let mut weyl = Array4::<f64>::zeros((4, 4, 4, 4));

//     for mu in 0..4 {
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     let mut term = riemann[[mu, nu, rho, sigma]];

//                     term -= 0.5 * (
//                         g[(mu, rho)] * ricci[(sigma, nu)] - g[(mu, sigma)] * ricci[(rho, nu)]
//                         - g[(nu, rho)] * ricci[(sigma, mu)] + g[(nu, sigma)] * ricci[(rho, mu)]
//                     );

//                     term += (scalar_r / 6.0)
//                         * (g[(mu, rho)] * g[(sigma, nu)] - g[(mu, sigma)] * g[(rho, nu)]);

//                     weyl[[mu, nu, rho, sigma]] = term;
//                 }
//             }
//         }
//     }
//     weyl
// }

// fn weyl_invariant(g_inv: &Matrix4<f64>, weyl: &Array4<f64>) -> f64 {
//     let mut sum = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     for alpha in 0..4 {
//                         for beta in 0..4 {
//                             for gamma in 0..4 {
//                                 for delta in 0..4 {
//                                     sum += g_inv[(mu, alpha)] * g_inv[(nu, beta)] *
//                                         g_inv[(rho, gamma)] * g_inv[(sigma, delta)] *
//                                         weyl[[mu, nu, rho, sigma]] * weyl[[alpha, beta, gamma, delta]];
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     sum
// }

// fn main() {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     println!("r\tC^2 (Weyl Invariant)");

//     for i in 1..100 {
//         let r = i as f64 * 0.1;

//         let g = metric(r, 0.0, r0, n, c);
//         let g_inv = g.try_inverse().unwrap();
//         let dg = metric_derivatives(r, r0, n, c, h);
//         let gamma = christoffel(&g, &g_inv, &dg);
//         let ricci = ricci_tensor(&gamma);
//         let scalar_r = scalar_curvature(&g_inv, &ricci);
//         let riemann = riemann_tensor(&gamma);
//         let weyl = weyl_tensor(&g, &g_inv, &ricci, scalar_r, &riemann);
//         let c2 = weyl_invariant(&g_inv, &weyl);

//         println!("{:.2}\t{:.6}", r, c2);
//     }
// }


// #11 тензор Бель-Робинсона

// use std::f64::consts::PI;
// use nalgebra::Matrix4;
// use ndarray::{Array3, Array4};

// fn metric(r: f64, _theta: f64, r0: f64, n: f64, c: f64) -> Matrix4<f64> {
//     let f = 1.0 / (1.0 + (r / r0).powf(n));
//     let g_rr = 16.0 * (1.0 + r).powi(6)
//         / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
//     let g_thetatheta = r.powi(2) / (1.0 + r).powi(2);

//     let mut g = Matrix4::<f64>::zeros();
//     g[(0, 0)] = -f * c * c;
//     g[(1, 1)] = g_rr;
//     g[(2, 2)] = g_thetatheta;
//     g[(3, 3)] = 1.0;
//     g
// }

// fn metric_derivatives(r: f64, r0: f64, n: f64, c: f64, h: f64) -> Matrix4<f64> {
//     let gp = metric(r + h, 0.0, r0, n, c);
//     let gm = metric(r - h, 0.0, r0, n, c);
//     let mut dg = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             dg[(mu, nu)] = (gp[(mu, nu)] - gm[(mu, nu)]) / (2.0 * h);
//         }
//     }
//     dg
// }

// fn christoffel(g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &Matrix4<f64>) -> Array3<f64> {
//     let mut gamma = Array3::<f64>::zeros((4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for nu in 0..4 {
//                 let mut sum = 0.0;
//                 for sigma in 0..4 {
//                     let dmu = if mu == 1 { dg[(nu, sigma)] } else { 0.0 };
//                     let dnu = if nu == 1 { dg[(mu, sigma)] } else { 0.0 };
//                     let dsigma = if sigma == 1 { dg[(mu, nu)] } else { 0.0 };
//                     sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
//                 }
//                 gamma[[lambda, mu, nu]] = 0.5 * sum;
//             }
//         }
//     }
//     gamma
// }

// fn riemann_tensor(gamma: &Array3<f64>) -> Array4<f64> {
//     let mut riemann = Array4::<f64>::zeros((4, 4, 4, 4));
//     for lambda in 0..4 {
//         for mu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     let mut term3 = 0.0;
//                     let mut term4 = 0.0;
//                     for nu in 0..4 {
//                         term3 += gamma[[lambda, nu, rho]] * gamma[[nu, mu, sigma]];
//                         term4 += gamma[[lambda, nu, sigma]] * gamma[[nu, mu, rho]];
//                     }
//                     riemann[[lambda, mu, rho, sigma]] = term3 - term4;
//                 }
//             }
//         }
//     }
//     riemann
// }

// fn ricci_tensor(gamma: &Array3<f64>) -> Matrix4<f64> {
//     let mut ricci = Matrix4::<f64>::zeros();
//     for mu in 0..4 {
//         for nu in 0..4 {
//             let mut sum = 0.0;
//             for lambda in 0..4 {
//                 sum += gamma[[lambda, mu, nu]] * gamma[[lambda, nu, mu]]
//                     - gamma[[lambda, lambda, mu]] * gamma[[nu, nu, lambda]];
//             }
//             ricci[(mu, nu)] = sum;
//         }
//     }
//     ricci
// }

// fn scalar_curvature(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
//     let mut r = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             r += g_inv[(mu, nu)] * ricci[(mu, nu)];
//         }
//     }
//     r
// }

// fn weyl_tensor(
//     g: &Matrix4<f64>,
//     g_inv: &Matrix4<f64>,
//     ricci: &Matrix4<f64>,
//     scalar_r: f64,
//     riemann: &Array4<f64>,
// ) -> Array4<f64> {
//     let mut weyl = Array4::<f64>::zeros((4, 4, 4, 4));

//     for mu in 0..4 {
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     let mut term = riemann[[mu, nu, rho, sigma]];

//                     term -= 0.5 * (
//                         g[(mu, rho)] * ricci[(sigma, nu)] - g[(mu, sigma)] * ricci[(rho, nu)]
//                         - g[(nu, rho)] * ricci[(sigma, mu)] + g[(nu, sigma)] * ricci[(rho, mu)]
//                     );

//                     term += (scalar_r / 6.0)
//                         * (g[(mu, rho)] * g[(sigma, nu)] - g[(mu, sigma)] * g[(rho, nu)]);

//                     weyl[[mu, nu, rho, sigma]] = term;
//                 }
//             }
//         }
//     }
//     weyl
// }

// fn weyl_invariant(g_inv: &Matrix4<f64>, weyl: &Array4<f64>) -> f64 {
//     let mut sum = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     for alpha in 0..4 {
//                         for beta in 0..4 {
//                             for gamma in 0..4 {
//                                 for delta in 0..4 {
//                                     sum += g_inv[(mu, alpha)] * g_inv[(nu, beta)] *
//                                         g_inv[(rho, gamma)] * g_inv[(sigma, delta)] *
//                                         weyl[[mu, nu, rho, sigma]] * weyl[[alpha, beta, gamma, delta]];
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     sum
// }

// fn belrobinson_tensor(
//     g: &Matrix4<f64>,
//     riemann: &Array4<f64>,
// ) -> f64 {
//     let mut belrobinson = 0.0;
//     for mu in 0..4 {
//         for nu in 0..4 {
//             for rho in 0..4 {
//                 for sigma in 0..4 {
//                     belrobinson += riemann[[mu, nu, rho, sigma]] * riemann[[mu, nu, rho, sigma]];
//                 }
//             }
//         }
//     }
//     belrobinson
// }

// fn main() {
//     let r0 = 1.0;
//     let n = 2.0;
//     let c = 1.0;
//     let h = 1e-5;

//     println!("r\tC^2 (Weyl Invariant)\tT_{{0000}} (Bel-Robinson)");

//     for i in 1..100 {
//         let r = i as f64 * 0.1;

//         let g = metric(r, 0.0, r0, n, c);
//         let g_inv = g.try_inverse().unwrap();
//         let dg = metric_derivatives(r, r0, n, c, h);
//         let gamma = christoffel(&g, &g_inv, &dg);
//         let ricci = ricci_tensor(&gamma);
//         let scalar_r = scalar_curvature(&g_inv, &ricci);
//         let riemann = riemann_tensor(&gamma);
//         let weyl = weyl_tensor(&g, &g_inv, &ricci, scalar_r, &riemann);
//         let c2 = weyl_invariant(&g_inv, &weyl);
//         let belrobinson = belrobinson_tensor(&g, &riemann);

//         println!("{:.2}\t{:.6}\t{:.6}", r, c2, belrobinson);
//     }
// }


//////// КВАНТОВЫЕ ТЕСТЫ

// use std::f64::consts::PI;

// // Основные физические константы
// const HBAR: f64 = 1.0545718e-34; // постоянная Планка / 2π, Дж·с
// const C: f64 = 299792458.0; // скорость света, м/с
// const G: f64 = 6.67430e-11; // гравитационная постоянная, м³/(кг·с²)

// // Параметры метрики HyperTwist
// const R0: f64 = 1.0;
// const N: f64 = 2.0;

// // Компоненты метрики
// fn f(r: f64) -> f64 {
//     1.0 / (1.0 + (r / R0).powf(N))
// }

// fn g_rr(r: f64) -> f64 {
//     16.0 * (1.0 + r).powi(6) / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2))
// }

// fn g_rtheta(r: f64) -> f64 {
//     -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4))
// }

// fn g_thetatheta(r: f64) -> f64 {
//     r.powi(2) / (1.0 + r).powi(2)
// }

// // Скалярная кривизна (аппроксимация на основе ваших данных)
// fn scalar_curvature(r: f64) -> f64 {
//     if r < 0.5 {
//         0.48 / r.powf(1.5) // аппроксимация для r < 0.5
//     } else {
//         0.13 / r.powi(2) // аппроксимация для r > 0.5
//     }
// }

// // Инвариант Вейля (аппроксимация на основе ваших данных)
// fn weyl_invariant(r: f64) -> f64 {
//     if r < 1.1 {
//         0.1 / r.powi(2) // положительный при r < 1.1
//     } else {
//         -0.1 / r.powi(2) // отрицательный при r > 1.1
//     }
// }

// // Квантовый потенциал
// fn quantum_potential(r: f64) -> f64 {
//     let r_curv = scalar_curvature(r);
//     // Для частицы массой 1e-27 кг
//     HBAR.powi(2) * r_curv / (2.0 * 1e-27 * C.powi(2))
// }

// // Плотность энергии вакуума
// fn vacuum_energy_density(r: f64) -> f64 {
//     let r_curv = scalar_curvature(r);
//     HBAR * C.powi(3) * r_curv / (8.0 * PI * G)
// }

// // Квантовая поправка к инварианту Вейля
// fn quantum_weyl_correction(r: f64) -> f64 {
//     let classical = weyl_invariant(r);
//     // Квантовая поправка в первом приближении (упрощенная модель)
//     let correction = HBAR * vacuum_energy_density(r) / (C.powi(3) * 1e20); // нормализация
//     classical + correction
// }

// // Расчет смещения точки изменения знака инварианта Вейля
// fn weyl_sign_change_shift() -> f64 {
//     // Теоретическая оценка смещения
//     HBAR * G / (C.powi(3) * R0.powi(2)) * scalar_curvature(1.1)
// }

// fn main() {
//     println!("# Анализ квантовых эффектов метрики HyperTwist");
//     println!("r\tf(r)\tg_rr(r)\tg_rtheta(r)\tg_thetatheta(r)");
    
//     let r_values = [0.1, 0.5, 1.0, 1.1, 2.0, 5.0, 10.0];
    
//     for &r in &r_values {
//         println!("{:.1}\t{:.6}\t{:.6}\t{:.6}\t{:.6}", 
//                  r, f(r), g_rr(r), g_rtheta(r), g_thetatheta(r));
//     }
    
//     println!("\n# Скалярная кривизна и квантовые эффекты");
//     println!("r\tR(r)\tV_квантовый\tρ_вакуума");
    
//     for &r in &r_values {
//         println!("{:.1}\t{:.6e}\t{:.6e}\t{:.6e}", 
//                  r, scalar_curvature(r), quantum_potential(r), vacuum_energy_density(r));
//     }
    
//     println!("\n# Инвариант Вейля и квантовые поправки");
//     println!("r\tC²_классич\tC²_квантовый\tИзменение(%)");
    
//     for &r in &r_values {
//         let c2_classic = weyl_invariant(r);
//         let c2_quantum = quantum_weyl_correction(r);
//         let change_pct = (c2_quantum - c2_classic) / c2_classic.abs() * 100.0;
        
//         println!("{:.1}\t{:.6e}\t{:.6e}\t{:.6}", 
//                  r, c2_classic, c2_quantum, change_pct);
//     }
    
//     println!("\n# Смещение точки изменения знака инварианта Вейля");
//     println!("Классическое значение: r = 1.1");
//     println!("Теоретическая оценка смещения: {:.6e}", weyl_sign_change_shift());
//     println!("С учетом квантовых поправок: r ≈ {:.6}", 1.1 + weyl_sign_change_shift());
    
//     // Оценка эффекта Казимира
//     println!("\n# Модификация силы Казимира в метрике HyperTwist");
//     println!("r\tF_плоское\tФактор_метрики\tF_модифиц");
    
//     for &r in &r_values {
//         let f_flat = -PI.powi(2) * HBAR * C / (240.0 * r.powi(4));
//         let metric_factor = (f(r) * g_rr(r) * g_thetatheta(r)).abs().sqrt();
//         let f_curved = f_flat * metric_factor;
        
//         println!("{:.1}\t{:.6e}\t{:.6}\t{:.6e}", r, f_flat, metric_factor, f_curved);
//     }
    
//     // Температура Хокинга для компактного объекта
//     println!("\n# Температура Хокинга для компактного объекта с метрикой HyperTwist");
//     println!("M/M_☉\tr_h\tT_стандарт(K)\tT_HyperTwist(K)\tОтношение");
    
//     let solar_mass = 1.989e30; // масса Солнца, кг
//     let masses = [1.0, 5.0, 10.0]; // массы в единицах солнечных масс
    
//     for &m_solar in &masses {
//         let mass = m_solar * solar_mass;
//         // Оценка радиуса горизонта (упрощенно)
//         let r_h = 2.0 * G * mass / C.powi(2);
        
//         // Стандартная температура Хокинга
//         let t_standard = HBAR * C.powi(3) / (8.0 * PI * G * mass * 1.380649e-23);
        
//         // Модифицированная температура с учетом скручивания
//         let alpha = 1.0; // коэффициент порядка единицы
//         let twist_factor = 1.0 + alpha * g_rtheta(r_h) / g_rr(r_h);
//         let t_hypertwist = t_standard * twist_factor;
        
//         println!("{:.1}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6}", 
//                  m_solar, r_h, t_standard, t_hypertwist, t_hypertwist / t_standard);
//     }
// }


/// #12 Доп тесты 
// use std::f64::consts::PI;
// use std::fs::File;
// use std::io::Write;

// // Космологические константы
// const H0: f64 = 70.0; // Постоянная Хаббла, км/с/Мпк
// const OM_M: f64 = 0.3; // Омега материи в стандартной модели
// const OM_L: f64 = 0.7; // Омега тёмной энергии в стандартной модели
// const C_LIGHT: f64 = 299792.458; // Скорость света, км/с

// // Параметры метрики HyperTwist
// const R0: f64 = 1.0;
// const N: f64 = 2.0;

// // Метрика HyperTwist
// fn f_metric(r: f64) -> f64 {
//     1.0 / (1.0 + (r / R0).powf(N))
// }

// fn g_rr(r: f64) -> f64 {
//     16.0 * (1.0 + r).powi(6) / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2))
// }

// fn g_rtheta(r: f64) -> f64 {
//     -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4))
// }

// fn g_thetatheta(r: f64) -> f64 {
//     r.powi(2) / (1.0 + r).powi(2)
// }

// // Космологическое расширение в модели ΛCDM
// fn hubble_lcdm(z: f64) -> f64 {
//     H0 * (OM_M * (1.0 + z).powi(3) + OM_L).sqrt()
// }

// fn comoving_distance_lcdm(z: f64, steps: usize) -> f64 {
//     // Интегрирование для получения сопутствующего расстояния
//     let dz = z / steps as f64;
//     let mut distance = 0.0;
    
//     for i in 0..steps {
//         let zi = i as f64 * dz;
//         distance += C_LIGHT / hubble_lcdm(zi) * dz;
//     }
    
//     distance
// }

// // Влияние скручивания на космологическое расширение (модель HyperTwist)
// fn hubble_hypertwist(z: f64, twist_effect: f64) -> f64 {
//     // twist_effect - параметр, определяющий влияние скручивания
//     // Для r < R0 (z > z0) скручивание усиливает расширение
//     // Для r > R0 (z < z0) скручивание замедляет расширение
    
//     let z0 = 1.1; // Соответствует r ≈ R0, где C² меняет знак
//     let twist_factor = if z < z0 {
//         1.0 - twist_effect * (z0 - z) / z0
//     } else {
//         1.0 + twist_effect * (z - z0) / (1.0 + z)
//     };
    
//     hubble_lcdm(z) * twist_factor
// }

// fn comoving_distance_hypertwist(z: f64, twist_effect: f64, steps: usize) -> f64 {
//     let dz = z / steps as f64;
//     let mut distance = 0.0;
    
//     for i in 0..steps {
//         let zi = i as f64 * dz;
//         distance += C_LIGHT / hubble_hypertwist(zi, twist_effect) * dz;
//     }
    
//     distance
// }

// // Эффективная плотность энергии скручивания как функция z
// fn twist_energy_density(z: f64) -> f64 {
//     let r = R0 / (1.0 + z); // Приближенное соответствие z и r
//     let twist_term = g_rtheta(r).powi(2) / (g_rr(r) * g_thetatheta(r));
    
//     // Нормализация для соответствия наблюдаемой темной энергии при z = 0
//     let norm = 0.7 / twist_term.abs();
//     twist_term * norm
// }

// // Расчет углового диаметра для стандартной свечи
// // Исправлено: теперь возвращает кортеж (f64, f64)
// fn angular_diameter(z: f64, physical_size: f64, twist_effect: f64) -> (f64, f64) {
//     let distance_lcdm = comoving_distance_lcdm(z, 1000) / (1.0 + z);
//     let distance_hypertwist = comoving_distance_hypertwist(z, twist_effect, 1000) / (1.0 + z);
    
//     let angle_lcdm = physical_size / distance_lcdm;
//     let angle_hypertwist = physical_size / distance_hypertwist;
    
//     (angle_hypertwist, angle_lcdm)
// }

// // Моделирование сигнала сверхновых типа Ia
// fn supernova_magnitude(z: f64, twist_effect: f64) -> (f64, f64) {
//     // Модуль расстояния в космологии
//     let distance_modulus_lcdm = 5.0 * (comoving_distance_lcdm(z, 1000) * (1.0 + z)).log10() + 25.0;
//     let distance_modulus_hypertwist = 5.0 * (comoving_distance_hypertwist(z, twist_effect, 1000) * (1.0 + z)).log10() + 25.0;
    
//     (distance_modulus_hypertwist, distance_modulus_lcdm)
// }

// // Возмущения в распределении материи
// fn matter_perturbation_growth(z: f64, twist_effect: f64) -> (f64, f64) {
//     // Упрощенная модель роста возмущений
//     let growth_lcdm = 1.0 / (1.0 + z);
    
//     // В метрике HyperTwist, скручивание усиливает гравитационное притяжение
//     let z0 = 1.1;
//     let twist_factor = if z < z0 {
//         1.0 + 0.5 * twist_effect * (z0 - z) / z0
//     } else {
//         1.0 - 0.3 * twist_effect * (z - z0) / (1.0 + z)
//     };
    
//     let growth_hypertwist = growth_lcdm * twist_factor;
    
//     (growth_hypertwist, growth_lcdm)
// }

// fn main() -> std::io::Result<()> {
//     let twist_effect = 0.3; // Параметр влияния скручивания (может требовать калибровки)
//     let z_values: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
    
//     // Файл для записи результатов Хаббла
//     let mut hubble_file = File::create("hubble_parameter.csv")?;
//     writeln!(hubble_file, "z,H_ΛCDM,H_HyperTwist,Twist_Energy_Density")?;
    
//     for &z in &z_values {
//         let h_lcdm = hubble_lcdm(z);
//         let h_hypertwist = hubble_hypertwist(z, twist_effect);
//         let twist_energy = twist_energy_density(z);
        
//         writeln!(hubble_file, "{},{},{},{}", z, h_lcdm, h_hypertwist, twist_energy)?;
//     }
    
//     // Файл для записи результатов сверхновых
//     let mut sn_file = File::create("supernova_magnitude.csv")?;
//     writeln!(sn_file, "z,m_ΛCDM,m_HyperTwist,Difference")?;
    
//     for &z in &z_values {
//         if z > 0.01 {  // Исключаем z ~ 0 для избежания особенностей
//             let (m_hypertwist, m_lcdm) = supernova_magnitude(z, twist_effect);
//             let diff = m_hypertwist - m_lcdm;
            
//             writeln!(sn_file, "{},{},{},{}", z, m_lcdm, m_hypertwist, diff)?;
//         }
//     }
    
//     // Файл для записи роста возмущений
//     let mut growth_file = File::create("perturbation_growth.csv")?;
//     writeln!(growth_file, "z,Growth_ΛCDM,Growth_HyperTwist,Ratio")?;
    
//     for &z in &z_values {
//         let (growth_hypertwist, growth_lcdm) = matter_perturbation_growth(z, twist_effect);
//         let ratio = growth_hypertwist / growth_lcdm;
        
//         writeln!(growth_file, "{},{},{},{}", z, growth_lcdm, growth_hypertwist, ratio)?;
//     }
    
//     // Файл для углового диаметра (Тест Алкока-Пачинского)
//     let mut angle_file = File::create("angular_diameter.csv")?;
//     writeln!(angle_file, "z,Angle_ΛCDM,Angle_HyperTwist,Ratio")?;
    
//     let physical_size = 1.0; // Размер в Мпк
//     for &z in &z_values {
//         if z > 0.01 {
//             let (angle_hypertwist, angle_lcdm) = angular_diameter(z, physical_size, twist_effect);
//             let ratio = angle_hypertwist / angle_lcdm;
            
//             writeln!(angle_file, "{},{},{},{}", z, angle_lcdm, angle_hypertwist, ratio)?;
//         }
//     }
    
//     println!("✅ Все расчеты завершены. Результаты записаны в CSV файлы.");
//     println!("📊 hubble_parameter.csv - параметр Хаббла и эффективная плотность энергии скручивания");
//     println!("📊 supernova_magnitude.csv - модули расстояния для сверхновых типа Ia");
//     println!("📊 perturbation_growth.csv - рост возмущений плотности материи");
//     println!("📊 angular_diameter.csv - тест углового диаметра (тест Алкока-Пачинского)");
    
//     Ok(())
// }


/// #13  - последние тесты
/// 
// use std::f64::consts::PI;
// use std::fs::File;
// use std::io::Write;
// use std::collections::HashMap;

// // Физические константы
// const C: f64 = 2.99792458e8; // скорость света, м/с
// const G: f64 = 6.67430e-11; // гравитационная постоянная, м³/(кг·с²)
// const MSUN: f64 = 1.989e30; // масса Солнца, кг

// // Параметры метрики HyperTwist
// const R0: f64 = 1.0;
// const N: f64 = 2.0;

// // Компоненты метрики
// fn f(r: f64) -> f64 {
//     1.0 / (1.0 + (r / R0).powf(N))
// }

// fn g_rr(r: f64) -> f64 {
//     16.0 * (1.0 + r).powi(6) / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2))
// }

// fn g_rtheta(r: f64) -> f64 {
//     -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4))
// }

// fn g_thetatheta(r: f64) -> f64 {
//     r.powi(2) / (1.0 + r).powi(2)
// }

// // Модифицированная скорость распространения гравитационных волн
// fn gw_speed(r: f64) -> f64 {
//     // В ОТО скорость гравитационных волн = скорости света
//     // В метрике HyperTwist скручивание может влиять на скорость распространения
    
//     let twist_factor = (g_rtheta(r) / g_rr(r)).powi(2) * r.powi(2);
    
//     // Эффект скручивания на скорость GW
//     // Нормализация для получения c при больших r
//     C * (1.0 + twist_factor / (1.0 + 10.0 * twist_factor))
// }

// // Дополнительная поляризационная мода из-за скручивания
// fn twist_polarization_amplitude(r: f64, base_amplitude: f64) -> f64 {
//     // В стандартной ОТО только два состояния поляризации (+, ×)
//     // В нашей метрике может появиться дополнительная мода из-за скручивания
    
//     let twist_factor = g_rtheta(r).abs() / (g_rr(r) * g_thetatheta(r)).sqrt();
    
//     // Амплитуда дополнительной поляризационной моды пропорциональна twist_factor
//     base_amplitude * twist_factor
// }

// // Вычисление искажения формы сигнала гравитационной волны от слияния черных дыр
// fn gw_waveform_distortion(r: f64, freq: f64) -> f64 {
//     // Вычисляет искажение фазы из-за дисперсии гравитационных волн
//     // в метрике HyperTwist относительно стандартной ОТО
    
//     // Частотная зависимость дисперсии из-за скручивания
//     let dispersion_factor = 2.0 * PI * freq * r / C;
    
//     // Искажение фазы как функция частоты и компоненты скручивания
//     dispersion_factor * g_rtheta(r).abs() / g_rr(r)
// }

// // Модель сливающихся черных дыр (упрощенная)
// struct BinaryMerger {
//     m1: f64,      // Масса первой черной дыры, солнечные массы
//     m2: f64,      // Масса второй черной дыры, солнечные массы
//     distance: f64, // Расстояние до наблюдателя, Мпк
// }

// impl BinaryMerger {
//     // Вычисляет базовую амплитуду гравитационной волны
//     fn base_amplitude(&self) -> f64 {
//         let mtot = (self.m1 + self.m2) * MSUN;
//         let dist_m = self.distance * 3.08567758e22; // Мпк в метры
        
//         // Формула базовой амплитуды GW в ОТО
//         (G * mtot) / (C * C * dist_m)
//     }
    
//     // Вычисляет характерную частоту слияния
//     fn merger_frequency(&self) -> f64 {
//         let mtot = (self.m1 + self.m2) * MSUN;
        
//         // Характерная частота в Гц
//         C.powi(3) / (G * mtot * PI * 6.0)
//     }
    
//     // Генерирует предсказанный LIGO сигнал в метриках ОТО и HyperTwist
//     fn generate_signal(&self, r_observer: f64) -> HashMap<String, Vec<(f64, f64)>> {
//         let base_amp = self.base_amplitude();
//         let f_merger = self.merger_frequency();
        
//         // Временной массив до и после слияния (в секундах)
//         let times: Vec<f64> = (-1000..1000).map(|i| i as f64 * 0.001).collect();
        
//         // Частотная эволюция (чирп)
//         let frequencies: Vec<f64> = times.iter()
//             .map(|&t| {
//                 if t < 0.0 {
//                     // До слияния - частота растет
//                     f_merger * (1.0 - t).powf(-0.375)
//                 } else {
//                     // После слияния - затухающий сигнал
//                     f_merger * (1.0 + t * 5.0).powf(-0.8)
//                 }
//             })
//             .collect();
        
//         // Амплитуды в ОТО
//         let amp_gtr: Vec<(f64, f64)> = times.iter().enumerate()
//             .map(|(i, &t)| {
//                 let phase = frequencies[i] * t * 2.0 * PI;
//                 let amp = if t < 0.0 {
//                     base_amp * (1.0 - t).powf(-0.25)
//                 } else {
//                     base_amp * (1.0 + t * 5.0).powf(-0.5)
//                 };
//                 (t, amp * phase.sin())
//             })
//             .collect();
        
//         // Амплитуды в HyperTwist
//         let amp_hypertwist: Vec<(f64, f64)> = times.iter().enumerate()
//             .map(|(i, &t)| {
//                 let freq = frequencies[i];
//                 let phase_distortion = gw_waveform_distortion(r_observer, freq);
//                 let phase = freq * t * 2.0 * PI + phase_distortion;
                
//                 let amp = if t < 0.0 {
//                     base_amp * (1.0 - t).powf(-0.25)
//                 } else {
//                     base_amp * (1.0 + t * 5.0).powf(-0.5)
//                 };
                
//                 // Дополнительная амплитуда из-за скручивания
//                 let twist_amp = twist_polarization_amplitude(r_observer, amp);
                
//                 // Суммарный сигнал из стандартной и дополнительной моды
//                 let total_amp = amp * phase.sin() + twist_amp * (phase * 2.0).sin();
                
//                 (t, total_amp)
//             })
//             .collect();
        
//         // Соберем результаты в словарь
//         let mut results = HashMap::new();
//         results.insert("GTR".to_string(), amp_gtr);
//         results.insert("HyperTwist".to_string(), amp_hypertwist);
//         results.insert("Frequencies".to_string(), 
//             frequencies.iter().enumerate().map(|(i, &f)| (times[i], f)).collect());
        
//         results
//     }
// }

// // Рассчитывает эффективное отношение сигнал/шум для двух моделей
// fn signal_to_noise_ratio(original: &[(f64, f64)], modified: &[(f64, f64)]) -> f64 {
//     let mut sum_diff_squared = 0.0;
//     let mut sum_orig_squared = 0.0;
    
//     for (i, &(_, amp_orig)) in original.iter().enumerate() {
//         let (_, amp_mod) = modified[i];
//         let diff = amp_orig - amp_mod;
        
//         sum_diff_squared += diff * diff;
//         sum_orig_squared += amp_orig * amp_orig;
//     }
    
//     if sum_orig_squared > 0.0 {
//         (sum_diff_squared / sum_orig_squared).sqrt()
//     } else {
//         0.0
//     }
// }

// // Анализ распространения гравитационных волн на разных расстояниях
// fn analyze_gw_propagation() -> std::io::Result<()> {
//     // Файл для скорости гравитационных волн
//     let mut speed_file = File::create("gw_speed.csv")?;
//     writeln!(speed_file, "r,GW_Speed_Standard,GW_Speed_HyperTwist,Ratio")?;
    
//     for i in 1..101 {
//         let r = i as f64 * 0.1;
//         let speed_ht = gw_speed(r);
//         let ratio = speed_ht / C;
        
//         writeln!(speed_file, "{},{},{},{}", r, C, speed_ht, ratio)?;
//     }
    
//     // Файл для дополнительных поляризационных мод
//     let mut pol_file = File::create("gw_polarization.csv")?;
//     writeln!(pol_file, "r,Standard_Modes,HyperTwist_Extra_Mode,Ratio")?;
    
//     let base_amp = 1.0e-21; // Типичная амплитуда GW
//     for i in 1..101 {
//         let r = i as f64 * 0.1;
//         let twist_amp = twist_polarization_amplitude(r, base_amp);
//         let ratio = twist_amp / base_amp;
        
//         writeln!(pol_file, "{},{},{},{}", r, base_amp, twist_amp, ratio)?;
//     }
    
//     Ok(())
// }

// // Сравнение детектируемости сигналов для разных типов двойных систем
// fn binary_merger_signals() -> std::io::Result<()> {
//     // Различные типы слияний черных дыр
//     let merger_types = vec![
//         BinaryMerger { m1: 10.0, m2: 10.0, distance: 400.0 },     // Равные массы, GW150914-подобное
//         BinaryMerger { m1: 30.0, m2: 30.0, distance: 1000.0 },    // Массивные ЧД
//         BinaryMerger { m1: 5.0, m2: 1.4, distance: 200.0 },       // ЧД + НЗ
//         BinaryMerger { m1: 1000.0, m2: 1000.0, distance: 3000.0 } // Сверхмассивные ЧД (LISA)
//     ];
    
//     let mut merger_file = File::create("binary_merger_snr.csv")?;
//     writeln!(merger_file, "BH1_Mass,BH2_Mass,Distance_Mpc,r_observer,GTR_Peak,HyperTwist_Peak,SNR_Diff")?;
    
//     for merger in &merger_types {
//         // Проверка на разных расстояниях от центра гравитирующего объекта
//         for r_obs in [0.5, 1.0, 2.0, 5.0, 10.0].iter() {
//             let signals = merger.generate_signal(*r_obs);
            
//             // Найдем максимальные амплитуды
//             let max_gtr = signals["GTR"].iter()
//                 .map(|&(_, amp)| amp.abs())
//                 .fold(0.0, f64::max);
            
//             let max_ht = signals["HyperTwist"].iter()
//                 .map(|&(_, amp)| amp.abs())
//                 .fold(0.0, f64::max);
            
//             // Вычислим эффективное отношение сигнал/шум различия моделей
//             let snr_diff = signal_to_noise_ratio(&signals["GTR"], &signals["HyperTwist"]);
            
//             writeln!(merger_file, "{},{},{},{},{:.3e},{:.3e},{:.6}",
//                      merger.m1, merger.m2, merger.distance, r_obs, max_gtr, max_ht, snr_diff)?;
//         }
//     }
    
//     // Подробное сравнение формы сигнала для выбранного слияния
//     let reference_merger = BinaryMerger { m1: 30.0, m2: 30.0, distance: 1000.0 };
//     let r_detailed = 1.0; // Радиус наблюдения вблизи топологического перехода
    
//     let signals = reference_merger.generate_signal(r_detailed);
    
//     let mut waveform_file = File::create("gw_waveform_comparison.csv")?;
//     writeln!(waveform_file, "Time,Frequency,GTR_Amplitude,HyperTwist_Amplitude,Difference")?;
    
//     for i in 0..signals["GTR"].len() {
//         let (t, amp_gtr) = signals["GTR"][i];
//         let (_, amp_ht) = signals["HyperTwist"][i];
//         let (_, freq) = signals["Frequencies"][i];
//         let diff = amp_ht - amp_gtr;
        
//         writeln!(waveform_file, "{},{},{},{},{}", t, freq, amp_gtr, amp_ht, diff)?;
//     }
    
//     println!("✅ Анализ гравитационных волн в метрике HyperTwist завершен.");
//     println!("📊 gw_speed.csv - Скорость распространения гравитационных волн");
//     println!("📊 gw_polarization.csv - Дополнительные поляризационные моды");
//     println!("📊 binary_merger_snr.csv - Сравнение сигналов слияния черных дыр");
//     println!("📊 gw_waveform_comparison.csv - Детальное сравнение формы сигнала");
    
//     Ok(())
// }

// fn main() -> std::io::Result<()> {
//     // Анализ распространения GW и поляризационных мод
//     analyze_gw_propagation()?;
    
//     // Сигналы от слияний двойных систем
//     binary_merger_signals()?;
    
//     Ok(())
// }



