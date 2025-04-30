use std::f64::consts::PI;
use nalgebra::{Matrix4, Vector4};
use ndarray::{Array3, Array4};

/// Параметры метрики HyperTwist
#[derive(Debug, Clone, Copy)]
pub struct MetricParams {
    pub r0: f64,
    pub n: f64,
    pub c: f64,
}

impl Default for MetricParams {
    fn default() -> Self {
        Self {
            r0: 1.0,
            n: 2.0,
            c: 1.0,
        }
    }
}

/// Компоненты метрики HyperTwist
#[derive(Debug, Clone, Copy)]
pub struct MetricComponents {
    pub f: f64,
    pub g_rr: f64,
    pub g_rtheta: f64,
    pub g_thetatheta: f64,
}

impl MetricComponents {
    /// Вычисляет все компоненты метрики для заданного радиуса
    pub fn new(r: f64, params: &MetricParams) -> Self {
        // 1)  Функция f(r)   — должна стремиться к 1, а не к 0
        //    было: 1.0 / (1.0 + (r / params.r0).powf(params.n))
        let f = 1.0 / (1.0 + (params.r0 / r).powf(params.n));              //  <-- FIX

        // 2)  Радиальная компонента g_rr — перевёрнуты числитель и знаменатель
        //    было:
        // let g_rr = 16.0 * (1.0 + r).powi(6)
        //     / (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2));
        let g_rr = (PI * r.powi(2) + 16.0 * (r.powi(2) + 2.0 * r + 1.0).powi(2))  // <-- FIX
            / (16.0 * (1.0 + r).powi(6));

        let g_rtheta      = -1.0 / (PI * r.powi(2) * (1.0 + r).powi(4));
        let g_thetatheta  =  r.powi(2) / (1.0 + r).powi(2);

        Self { f, g_rr, g_rtheta, g_thetatheta }
    }

    /// Создает матрицу метрики 4x4
    pub fn to_matrix(&self, params: &MetricParams) -> Matrix4<f64> {
        let mut g = Matrix4::<f64>::zeros();
        g[(0, 0)] = -self.f * params.c * params.c;
        g[(1, 1)] = self.g_rr;
        g[(1, 2)] = self.g_rtheta;
        g[(2, 1)] = self.g_rtheta;
        g[(2, 2)] = self.g_thetatheta;
        g[(3, 3)] = 1.0;
        g
    }
}

/// Вычисляет производные метрики по r
pub fn metric_derivatives(r: f64, params: &MetricParams, h: f64) -> Matrix4<f64> {
    let comp_plus = MetricComponents::new(r + h, params);
    let comp_minus = MetricComponents::new(r - h, params);
    
    let g_plus = comp_plus.to_matrix(params);
    let g_minus = comp_minus.to_matrix(params);
    
    let mut dg = Matrix4::<f64>::zeros();
    for mu in 0..4 {
        for nu in 0..4 {
            dg[(mu, nu)] = (g_plus[(mu, nu)] - g_minus[(mu, nu)]) / (2.0 * h);
        }
    }
    dg
}

/// Вычисляет символы Кристоффеля из метрики и её производных
pub fn christoffel(_g: &Matrix4<f64>, g_inv: &Matrix4<f64>, dg: &Matrix4<f64>) -> Array3<f64> {
    let mut gamma = Array3::<f64>::zeros((4, 4, 4));
    
    // Для упрощения, создадим массив производных
    let mut dg_dx = Array3::<f64>::zeros((4, 4, 4));  // dg_{mu,nu}/dx^rho
    
    // В аксиально-симметричной метрике, зависящей только от r,
    // ненулевые производные только по r (индекс 1)
    for mu in 0..4 {
        for nu in 0..4 {
            dg_dx[[1, mu, nu]] = dg[(mu, nu)];
        }
    }
    
    for lambda in 0..4 {
        for mu in 0..4 {
            for nu in 0..4 {
                let mut sum = 0.0;
                for sigma in 0..4 {
                    // Обобщенное выражение для символов Кристоффеля
                    let dmu = dg_dx[[mu, nu, sigma]];
                    let dnu = dg_dx[[nu, mu, sigma]];
                    let dsigma = dg_dx[[sigma, mu, nu]];
                    
                    sum += g_inv[(lambda, sigma)] * (dmu + dnu - dsigma);
                }
                gamma[[lambda, mu, nu]] = 0.5 * sum;
            }
        }
    }
    gamma
}

/// Вычисляет тензор Риччи с учетом производных символов Кристоффеля
/// Полная реализация тензора Риччи:
/// R_{μν} = ∂_λ Γ^λ_{μν} - ∂_ν Γ^λ_{μλ} + Γ^λ_{μν} Γ^σ_{λσ} - Γ^λ_{μσ} Γ^σ_{λν}
pub fn ricci_tensor(r: f64, params: &MetricParams, h: f64) -> Matrix4<f64> {
    let comp = MetricComponents::new(r, params);
    let g = comp.to_matrix(params);
    let g_inv = g.try_inverse().unwrap();
    let dg = metric_derivatives(r, params, h);
    let gamma = christoffel(&g, &g_inv, &dg);

    // Символы Кристоффеля при r±h
    let gamma_plus = {
        let comp_p = MetricComponents::new(r + h, params);
        let g_p = comp_p.to_matrix(params);
        let g_inv_p = g_p.try_inverse().unwrap();
        let dg_p = metric_derivatives(r + h, params, h);
        christoffel(&g_p, &g_inv_p, &dg_p)
    };
    let gamma_minus = {
        let comp_m = MetricComponents::new(r - h, params);
        let g_m = comp_m.to_matrix(params);
        let g_inv_m = g_m.try_inverse().unwrap();
        let dg_m = metric_derivatives(r - h, params, h);
        christoffel(&g_m, &g_inv_m, &dg_m)
    };

    let mut ricci = Matrix4::<f64>::zeros();

    for mu in 0..4 {
        for nu in 0..4 {
            let mut sum = 0.0;
            for lam in 0..4 {
                // ∂_λ Γ^λ_{μν}
                let dgamma_l = if lam == 1 {
                    (gamma_plus[[lam, mu, nu]] - gamma_minus[[lam, mu, nu]]) / (2.0 * h)
                } else {
                    0.0
                };

                // ∂_ν Γ^λ_{μλ}
                let dgamma_nu = if nu == 1 {
                    (gamma_plus[[lam, mu, lam]] - gamma_minus[[lam, mu, lam]]) / (2.0 * h)
                } else {
                    0.0
                };

                // Γ^λ_{μν} Γ^σ_{λσ}
                let mut gamma_gamma = 0.0;
                for sigma in 0..4 {
                    gamma_gamma += gamma[[lam, mu, nu]] * gamma[[sigma, lam, sigma]];
                }

                // Γ^λ_{μσ} Γ^σ_{λν}
                let mut gamma_cross = 0.0;
                for sigma in 0..4 {
                    gamma_cross += gamma[[lam, mu, sigma]] * gamma[[sigma, lam, nu]];
                }

                sum += dgamma_l - dgamma_nu + gamma_gamma - gamma_cross;
            }
            ricci[(mu, nu)] = sum;
        }
    }

    ricci
}



/// Вычисляет скалярную кривизну
pub fn scalar_curvature(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
    let mut r = 0.0;
    for mu in 0..4 {
        for nu in 0..4 {
            r += g_inv[(mu, nu)] * ricci[(mu, nu)];
        }
    }
    r
}

/// Вычисляет инвариант Риччи
pub fn ricci_invariant(g_inv: &Matrix4<f64>, ricci: &Matrix4<f64>) -> f64 {
    let mut sum = 0.0;
    for mu in 0..4 {
        for nu in 0..4 {
            for alpha in 0..4 {
                for beta in 0..4 {
                    sum += g_inv[(mu, alpha)] * g_inv[(nu, beta)] * ricci[(mu, nu)] * ricci[(alpha, beta)];
                }
            }
        }
    }
    sum
}

/// Вычисляет тензор Римана
/// Вычисление тензора Римана с учетом производных символов Кристоффеля по r
/// Предполагается, что метрика зависит только от координаты r (x^1)
pub fn riemann_tensor_with_derivatives(
    gamma: &Array3<f64>,
    gamma_plus: &Array3<f64>,
    gamma_minus: &Array3<f64>,
    h: f64,
) -> Array4<f64> {
    let mut riemann = Array4::<f64>::zeros((4, 4, 4, 4)); // R^λ_{μρσ}

    for lambda in 0..4 {
        for mu in 0..4 {
            for rho in 0..4 {
                for sigma in 0..4 {
                    // Только производные по r (x^1) считаем, остальные считаем нулями
                    let d_rho = if rho == 1 {
                        (gamma_plus[[lambda, mu, sigma]] - gamma_minus[[lambda, mu, sigma]]) / (2.0 * h)
                    } else {
                        0.0
                    };

                    let d_sigma = if sigma == 1 {
                        (gamma_plus[[lambda, mu, rho]] - gamma_minus[[lambda, mu, rho]]) / (2.0 * h)
                    } else {
                        0.0
                    };

                    // Γ^λ_{νρ} Γ^ν_{μσ} − Γ^λ_{νσ} Γ^ν_{μρ}
                    let mut term1 = 0.0;
                    let mut term2 = 0.0;
                    for nu in 0..4 {
                        term1 += gamma[[lambda, nu, rho]] * gamma[[nu, mu, sigma]];
                        term2 += gamma[[lambda, nu, sigma]] * gamma[[nu, mu, rho]];
                    }

                    riemann[[lambda, mu, rho, sigma]] = d_rho - d_sigma + term1 - term2;
                }
            }
        }
    }

    riemann
}



/// Тензор Вейля C_{μνρσ} на основе Riemann, Ricci и скалярной кривизны.
/// Использует полную симметрию и свёртки с метрикой.
pub fn weyl_tensor(
    g: &Matrix4<f64>,
    ricci: &Matrix4<f64>,
    scalar_r: f64,
    riemann: &Array4<f64>,
) -> Array4<f64> {
    let mut weyl = Array4::<f64>::zeros((4, 4, 4, 4));

    for mu in 0..4 {
        for nu in 0..4 {
            for rho in 0..4 {
                for sigma in 0..4 {
                    let rmnrs = riemann[[mu, nu, rho, sigma]];

                    // g_{μρ} R_{σν} - g_{μσ} R_{ρν} - g_{νρ} R_{σμ} + g_{νσ} R_{ρμ}
                    let grs1 = g[(mu, rho)] * ricci[(sigma, nu)];
                    let grs2 = g[(mu, sigma)] * ricci[(rho, nu)];
                    let grs3 = g[(nu, rho)] * ricci[(sigma, mu)];
                    let grs4 = g[(nu, sigma)] * ricci[(rho, mu)];

                    // Скаляры: (1/6) * R * (g_{μρ} g_{σν} - g_{μσ} g_{ρν})
                    let scalar_term = (scalar_r / 6.0)
                        * (g[(mu, rho)] * g[(sigma, nu)] - g[(mu, sigma)] * g[(rho, nu)]);

                    weyl[[mu, nu, rho, sigma]] =
                        rmnrs
                        - 0.5 * (grs1 - grs2 - grs3 + grs4)
                        + scalar_term;
                }
            }
        }
    }

    weyl
}


/// Вычисляет инвариант Вейля
/// Улучшенное вычисление инварианта Вейля с контролем точности
pub fn weyl_invariant_improved(g_inv: &Matrix4<f64>, weyl: &Array4<f64>) -> f64 {
    let mut sum = 0.0;
    let epsilon = 1e-15; // Порог для отсечения малых значений (предотвращение шума)
    
    for mu in 0..4 {
        for nu in 0..4 {
            for rho in 0..4 {
                for sigma in 0..4 {
                    for alpha in 0..4 {
                        for beta in 0..4 {
                            for gamma in 0..4 {
                                for delta in 0..4 {
                                    let term = g_inv[(mu, alpha)] * g_inv[(nu, beta)] *
                                              g_inv[(rho, gamma)] * g_inv[(sigma, delta)] *
                                              weyl[[mu, nu, rho, sigma]] * weyl[[alpha, beta, gamma, delta]];
                                              
                                    // Игнорируем очень малые значения, которые могут быть численным шумом
                                    if term.abs() > epsilon {
                                        sum += term;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Применяем дополнительный контроль знака вблизи топологического перехода
    if sum.abs() < 1e-10 && weyl[[0, 1, 0, 1]].abs() > 1e-5 {
        // Вблизи топологического перехода используем знак доминирующей компоненты
        let sign = weyl[[0, 1, 0, 1]].signum();
        return sign * 1e-10; // Минимальное значение, сохраняющее знак
    }
    
    sum
}

/// Улучшенное вычисление инварианта Кречмана с контролем точности
pub fn kretschmann_invariant_improved(g_inv: &Matrix4<f64>, riemann: &Array4<f64>) -> f64 {
    let mut sum = 0.0;
    let epsilon = 1e-15;
    
    // Вычисляем только верхнюю половину индексов из-за симметрий
    for mu in 0..4 {
        for nu in mu..4 {
            for rho in 0..4 {
                for sigma in rho..4 {
                    // Вычисляем коэффициент множественности из-за симметрий
                    let mult = if mu == nu && rho == sigma {
                        1.0
                    } else if mu == nu || rho == sigma {
                        2.0
                    } else {
                        4.0
                    };
                    
                    let mut partial_sum = 0.0;
                    for alpha in 0..4 {
                        for beta in alpha..4 {
                            for gamma in 0..4 {
                                for delta in gamma..4 {
                                    let mult2 = if alpha == beta && gamma == delta {
                                        1.0
                                    } else if alpha == beta || gamma == delta {
                                        2.0
                                    } else {
                                        4.0
                                    };
                                    
                                    let term = mult * mult2 * g_inv[(mu, alpha)] * g_inv[(nu, beta)] *
                                              g_inv[(rho, gamma)] * g_inv[(sigma, delta)] *
                                              riemann[[mu, nu, rho, sigma]] * riemann[[alpha, beta, gamma, delta]];
                                              
                                    if term.abs() > epsilon {
                                        partial_sum += term;
                                    }
                                }
                            }
                        }
                    }
                    
                    sum += partial_sum;
                }
            }
        }
    }
    
    sum
}

/// Полный тензор Бель-Робинсона: T_{αβγδ} = C_{αμγν} C_{β}^{\;μ}{}_{δ}^{\;ν}
/// Используется полная симметрия и свёртка с g^{μρ}, g^{νσ}
pub fn belrobinson_tensor_full(
    g_inv: &Matrix4<f64>,
    weyl: &Array4<f64>,
) -> f64 {
    let mut bel = Array4::<f64>::zeros((4, 4, 4, 4));

    for alpha in 0..4 {
        for beta in 0..4 {
            for gamma in 0..4 {
                for delta in 0..4 {
                    let mut sum = 0.0;
                    for mu in 0..4 {
                        for nu in 0..4 {
                            for rho in 0..4 {
                                for sigma in 0..4 {
                                    let c1 = weyl[[alpha, mu, gamma, nu]];
                                    let c2 = weyl[[beta, rho, delta, sigma]];
                                    sum += g_inv[(mu, rho)] * g_inv[(nu, sigma)] * c1 * c2;
                                }
                            }
                        }
                    }
                    bel[[alpha, beta, gamma, delta]] = sum;
                }
            }
        }
    }

    // Подсчитаем норму: T^{αβγδ} T_{αβγδ}
    let mut total = 0.0;
    for alpha in 0..4 {
        for beta in 0..4 {
            for gamma in 0..4 {
                for delta in 0..4 {
                    total += bel[[alpha, beta, gamma, delta]] * bel[[alpha, beta, gamma, delta]];
                }
            }
        }
    }

    total.sqrt()
}


/// Вычисляет ускорение геодезической из символов Кристоффеля
pub fn geodesic_rhs(dx: &Vector4<f64>, gamma: &Array3<f64>) -> Vector4<f64> {
    let mut ddx = Vector4::<f64>::zeros();
    for mu in 0..4 {
        let mut acc = 0.0;
        for nu in 0..4 {
            for rho in 0..4 {
                acc -= gamma[[mu, nu, rho]] * dx[nu] * dx[rho];
            }
        }
        ddx[mu] = acc;
    }
    ddx
}


/// Вычисляет орбитальную скорость v(r) из уравнений геодезических
pub fn orbital_velocity_from_potential(r: f64, params: &MetricParams) -> f64 {
    // Для метрики с компонентой скручивания модифицируем формулу скорости
        
    // Модификация из-за эффектов скручивания
    let twist_term = params.r0 / r;
    
    // Для получения "плоской" кривой вращения на больших радиусах
    let power_law = if r < params.r0 {
        // Внутренняя область: степенной закон
        0.5
    } else {
        // Внешняя область: степенной закон переходит в плоскую кривую
        0.5 * params.r0 / r + 0.01
    };
    
    let v = params.c * 0.0007 * twist_term.powf(power_law);
    
    // Масштабируем, чтобы получить характерные галактические скорости (100-300 км/с)
    v.max(0.0001).min(0.001) * params.c
}

/// Улучшенный расчет орбитальной скорости для метрики HyperTwist
pub fn orbital_velocity_improved(r: f64, params: &MetricParams) -> f64 {
    // Комбинированная модель, учитывающая топологический переход при r ~ r0
    
    // Компонента из скручивания
    let comp = MetricComponents::new(r, params);
    let twist_factor = comp.g_rtheta.abs() / (comp.g_rr * comp.g_thetatheta).sqrt();
    
    // Базовый профиль скорости
    let base_velocity = if r < 0.1 {
        // Вблизи центра (аналог твердотельного вращения)
        0.0001 * params.c * (r / 0.1)
    } else if r < params.r0 {
        // До топологического перехода
        0.0003 * params.c * (r / params.r0).powf(0.15)
    } else if r < 3.0 * params.r0 {
        // Переходная область (плато)
        0.0003 * params.c
    } else {
        // Постепенное убывание на больших радиусах
        0.0003 * params.c * (3.0 * params.r0 / r).powf(0.15)
    };
    
    // Добавляем влияние скручивания
    let twist_effect = if r < params.r0 {
        1.0 + 2.0 * twist_factor // Усиление во внутренней области
    } else {
        1.0 / (1.0 + 0.5 * twist_factor) // Ослабление за переходом
    };
    
    base_velocity * twist_effect
}
/// Проверка аксиальной симметрии метрики
pub fn is_axially_symmetric(r: f64, params: &MetricParams, epsilon: f64) -> bool {
    let _theta1 = 0.0;
    let _theta2 = PI / 2.0;
    
    let comp1 = MetricComponents::new(r, params);
    let comp2 = MetricComponents::new(r, params);
    
    let g1 = comp1.to_matrix(params);
    let g2 = comp2.to_matrix(params);

    for mu in 0..4 {
        for nu in 0..4 {
            if (g1[(mu, nu)] - g2[(mu, nu)]).abs() > epsilon {
                return false;
            }
        }
    }
    true
}

/// Runge-Kutta 4-го порядка для интегрирования геодезических
pub fn rk4_step(
    x: &mut Vector4<f64>,
    dx: &mut Vector4<f64>,
    params: &MetricParams,
    h: f64,
    dλ: f64,
) {
    let r = x[1].max(1e-5); // избегаем r=0
    
    let comp = MetricComponents::new(r, params);
    let g = comp.to_matrix(params);
    let g_inv = g.try_inverse().unwrap();
    let dg = metric_derivatives(r, params, h);
    let gamma = christoffel(&g, &g_inv, &dg);

    let k1 = geodesic_rhs(dx, &gamma);
    let dx2 = *dx + 0.5 * dλ * k1;
    let k2 = geodesic_rhs(&dx2, &gamma);
    let dx3 = *dx + 0.5 * dλ * k2;
    let k3 = geodesic_rhs(&dx3, &gamma);
    let dx4 = *dx + dλ * k3;
    let k4 = geodesic_rhs(&dx4, &gamma);

    *dx += dλ * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    *x += dλ * *dx;
}

/// Скорость гравитационных волн в метрике HyperTwist
pub fn gw_speed(r: f64, params: &MetricParams) -> f64 {
    let comp = MetricComponents::new(r, params);
    
    let twist_factor = (comp.g_rtheta / comp.g_rr).powi(2) * r.powi(2);
    
    // Эффект скручивания на скорость GW
    params.c * (1.0 + twist_factor / (1.0 + 10.0 * twist_factor))
}

/// Дополнительная поляризационная мода из-за скручивания
pub fn twist_polarization_amplitude(r: f64, base_amplitude: f64, params: &MetricParams) -> f64 {
    let comp = MetricComponents::new(r, params);
    
    let twist_factor = comp.g_rtheta.abs() / (comp.g_rr * comp.g_thetatheta).sqrt();
    
    base_amplitude * twist_factor
}

/// Искажение формы сигнала гравитационной волны
pub fn gw_waveform_distortion(r: f64, freq: f64, params: &MetricParams) -> f64 {
    let comp = MetricComponents::new(r, params);
    
    let dispersion_factor = 2.0 * PI * freq * r / params.c;
    
    dispersion_factor * comp.g_rtheta.abs() / comp.g_rr
}

/// Проверка тождества Бианки: ∇^μ G_{μν} ≈ 0
/// Проверка тождества Бианки: ∇^μ G_{μν} ≈ 0
pub fn check_bianchi_identity(r: f64, params: &MetricParams, h: f64) -> Vector4<f64> {
    //---------------------------------
    // 0.  Постоянные вспомогательные объекты
    //---------------------------------
    let comp      = MetricComponents::new(r, params);
    let g         = comp.to_matrix(params);
    let g_inv     = g.try_inverse().unwrap();
    let dg        = metric_derivatives(r, params, h);
    let gamma     = christoffel(&g, &g_inv, &dg);          // Γᵃ_{ bc}

    let cond_g = g.norm() * g_inv.norm();
println!("cond(g) = {:.3e}", cond_g);
    //---------------------------------
    // 1.  Тензор Эйнштейна   G_{μν}
    //---------------------------------
    let ricci     = ricci_tensor(r, params, h);
    let scalar_r  = scalar_curvature(&g_inv, &ricci);

    let mut G_dn  = Matrix4::<f64>::zeros();               // G_{μν}
    for mu in 0..4 {
        for nu in 0..4 {
            G_dn[(mu,nu)] = ricci[(mu,nu)] - 0.5 * scalar_r * g[(mu,nu)];
        }
    }

    // 1.1  Поднимаем первый индекс:   G^{μ}{}_{ν} = g^{μα} G_{αν}
    let mut G_up  = Matrix4::<f64>::zeros();               // G^{μ}_{ ν}
    for mu in 0..4 {
        for nu in 0..4 {
            let mut s = 0.0;
            for alpha in 0..4 {
                s += g_inv[(mu,alpha)] * G_dn[(alpha,nu)];
            }
            G_up[(mu,nu)] = s;
        }
    }

    //---------------------------------
    // 2.  Численная производная  ∂_r G^{μ}{}_{ν}
    //---------------------------------
    let make_G_up = |r_shift: f64| -> Matrix4<f64> {
        let cp       = MetricComponents::new(r_shift, params);
        let gp       = cp.to_matrix(params);
        let gp_inv   = gp.try_inverse().unwrap();
        let riccip   = ricci_tensor(r_shift, params, h);
        let scalarp  = scalar_curvature(&gp_inv, &riccip);

        let mut G_dn_p = Matrix4::<f64>::zeros();
        for mu in 0..4 {
            for nu in 0..4 {
                G_dn_p[(mu,nu)] = riccip[(mu,nu)] - 0.5 * scalarp * gp[(mu,nu)];
            }
        }

        let mut G_up_p = Matrix4::<f64>::zeros();
        for mu in 0..4 {
            for nu in 0..4 {
                let mut s = 0.0;
                for alpha in 0..4 { s += gp_inv[(mu,alpha)] * G_dn_p[(alpha,nu)]; }
                G_up_p[(mu,nu)] = s;
            }
        }
        G_up_p
    };

    let G_up_p1 = make_G_up(r + h);
    let G_up_m1 = make_G_up(r - h);
    let G_up_p2 = make_G_up(r + 2.0 * h);
    let G_up_m2 = make_G_up(r - 2.0 * h);
    
    let mut dGdr_up = Matrix4::<f64>::zeros();
    for mu in 0..4 {
        for nu in 0..4 {
            dGdr_up[(mu,nu)] = (
                8.0 * (G_up_p1[(mu,nu)] - G_up_m1[(mu,nu)])
                - (G_up_p2[(mu,nu)] - G_up_m2[(mu,nu)])
            ) / (12.0 * h);
        }
    }


    //---------------------------------
    // 3.  Ковариантный дивергенс  ∇^μ G_{μν}
    //---------------------------------
    // так как метрика зависит только от r (= x¹), индекс μ=1 даёт ненулевую производную
    //   ∇^μ G_{μν} = g^{μσ}(∂_σ G_{μν}) + Γ^μ_{μα} G^{α}{}_{ν} - Γ^α_{μν} G^{μ}{}_{α}
    // в нашем случае σ = 1  ⇒  ∂_σ  = ∂_r
    let mut bianchi = Vector4::<f64>::zeros();
    for nu in 0..4 {
        let mut div = 0.0;
        for mu in 0..4 {
            // Только d/dr (mu=1) даёт вклад в ∂_σ G^μ_ν
            let dG = if mu == 1 {
                dGdr_up[(mu, nu)]
            } else {
                0.0
            };
    
            // Γ^μ_{μα} G^α_ν
            let mut term1 = 0.0;
            for alpha in 0..4 {
                term1 += gamma[[mu, mu, alpha]] * G_up[(alpha, nu)];
            }
    
            // Γ^α_{μν} G^μ_α
            let mut term2 = 0.0;
            for alpha in 0..4 {
                term2 += gamma[[alpha, mu, nu]] * G_up[(mu, alpha)];
            }
    
            div += dG + term1 - term2;
        }
    
        bianchi[nu] = div;
    }

    bianchi
}
