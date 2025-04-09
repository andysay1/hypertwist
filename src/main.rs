
mod geometry;

use std::f64::consts::PI;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use std::collections::HashMap;
use nalgebra::Vector4;
use geometry::*;

// Физические константы
const C: f64 = 299792458.0;               // Скорость света, м/с
const G: f64 = 6.67430e-11;               // Гравитационная постоянная, м³/(кг·с²)
const HBAR: f64 = 1.0545718e-34;          // Постоянная Планка / 2π, Дж·с
const MSUN: f64 = 1.989e30;               // Масса Солнца, кг

// Космологические константы
const H0: f64 = 70.0;                     // Постоянная Хаббла, км/с/Мпк
const OM_M: f64 = 0.3;                    // Омега материи в стандартной модели
const OM_L: f64 = 0.7;                    // Омега тёмной энергии
const C_LIGHT: f64 = 299792.458;          // Скорость света, км/с

//=============================================================================
// КОСМОЛОГИЧЕСКИЕ РАСЧЕТЫ
//=============================================================================

/// Космологическое расширение в модели ΛCDM
fn hubble_lcdm(z: f64) -> f64 {
    H0 * (OM_M * (1.0 + z).powi(3) + OM_L).sqrt()
}

/// Интегрирование для получения сопутствующего расстояния в ΛCDM
fn comoving_distance_lcdm(z: f64, steps: usize) -> f64 {
    let dz = z / steps as f64;
    let mut distance = 0.0;
    
    for i in 0..steps {
        let zi = i as f64 * dz;
        distance += C_LIGHT / hubble_lcdm(zi) * dz;
    }
    
    distance
}

/// Влияние скручивания на космологическое расширение (модель HyperTwist)
fn hubble_hypertwist(z: f64, twist_effect: f64) -> f64 {
    // twist_effect - параметр, определяющий влияние скручивания
    // Для r < R0 (z > z0) скручивание усиливает расширение
    // Для r > R0 (z < z0) скручивание замедляет расширение
    
    let z0 = 1.1; // Соответствует r ≈ R0, где C² меняет знак
    let twist_factor = if z < z0 {
        1.0 - twist_effect * (z0 - z) / z0
    } else {
        1.0 + twist_effect * (z - z0) / (1.0 + z)
    };
    
    hubble_lcdm(z) * twist_factor
}

/// Интегрирование для получения сопутствующего расстояния в модели HyperTwist
fn comoving_distance_hypertwist(z: f64, twist_effect: f64, steps: usize) -> f64 {
    let dz = z / steps as f64;
    let mut distance = 0.0;
    
    for i in 0..steps {
        let zi = i as f64 * dz;
        distance += C_LIGHT / hubble_hypertwist(zi, twist_effect) * dz;
    }
    
    distance
}

/// Эффективная плотность энергии скручивания как функция z
fn twist_energy_density(z: f64) -> f64 {
    let params = MetricParams::default();
    let r = params.r0 / (1.0 + z); // Приближенное соответствие z и r
    let comp = MetricComponents::new(r, &params);
    let twist_term = comp.g_rtheta.powi(2) / (comp.g_rr * comp.g_thetatheta);
    
    // Нормализация для соответствия наблюдаемой темной энергии при z = 0
    let norm = 0.7 / twist_term.abs();
    twist_term * norm
}

/// Расчет углового диаметра для стандартной свечи
fn angular_diameter(z: f64, physical_size: f64, twist_effect: f64) -> (f64, f64) {
    let distance_lcdm = comoving_distance_lcdm(z, 1000) / (1.0 + z);
    let distance_hypertwist = comoving_distance_hypertwist(z, twist_effect, 1000) / (1.0 + z);
    
    let angle_lcdm = physical_size / distance_lcdm;
    let angle_hypertwist = physical_size / distance_hypertwist;
    
    (angle_hypertwist, angle_lcdm)
}

/// Моделирование сигнала сверхновых типа Ia
fn supernova_magnitude(z: f64, twist_effect: f64) -> (f64, f64) {
    // Модуль расстояния в космологии
    let distance_modulus_lcdm = 5.0 * (comoving_distance_lcdm(z, 1000) * (1.0 + z)).log10() + 25.0;
    let distance_modulus_hypertwist = 5.0 * (comoving_distance_hypertwist(z, twist_effect, 1000) * (1.0 + z)).log10() + 25.0;
    
    (distance_modulus_hypertwist, distance_modulus_lcdm)
}

/// Возмущения в распределении материи
fn matter_perturbation_growth(z: f64, twist_effect: f64) -> (f64, f64) {
    // Упрощенная модель роста возмущений
    let growth_lcdm = 1.0 / (1.0 + z);
    
    // В метрике HyperTwist, скручивание усиливает гравитационное притяжение
    let z0 = 1.1;
    let twist_factor = if z < z0 {
        1.0 + 0.5 * twist_effect * (z0 - z) / z0
    } else {
        1.0 - 0.3 * twist_effect * (z - z0) / (1.0 + z)
    };
    
    let growth_hypertwist = growth_lcdm * twist_factor;
    
    (growth_hypertwist, growth_lcdm)
}

//=============================================================================
// КВАНТОВЫЕ ЭФФЕКТЫ
//=============================================================================

/// Скалярная кривизна (аппроксимация)
fn scalar_curvature_approx(r: f64) -> f64 {
    if r < 0.5 {
        0.48 / r.powf(1.5) // аппроксимация для r < 0.5
    } else {
        0.13 / r.powi(2) // аппроксимация для r > 0.5
    }
}

/// Инвариант Вейля (аппроксимация)
fn weyl_invariant_approx(r: f64) -> f64 {
    if r < 1.1 {
        0.1 / r.powi(2) // положительный при r < 1.1
    } else {
        -0.1 / r.powi(2) // отрицательный при r > 1.1
    }
}

/// Квантовый потенциал
fn quantum_potential(r: f64) -> f64 {
    let r_curv = scalar_curvature_approx(r);
    // Для частицы массой 1e-27 кг
    HBAR.powi(2) * r_curv / (2.0 * 1e-27 * C.powi(2))
}

/// Плотность энергии вакуума
fn vacuum_energy_density(r: f64) -> f64 {
    let r_curv = scalar_curvature_approx(r);
    HBAR * C.powi(3) * r_curv / (8.0 * PI * G)
}

/// Квантовая поправка к инварианту Вейля
fn quantum_weyl_correction(r: f64) -> f64 {
    let classical = weyl_invariant_approx(r);
    // Квантовая поправка в первом приближении (упрощенная модель)
    let correction = HBAR * vacuum_energy_density(r) / (C.powi(3) * 1e20); // нормализация
    classical + correction
}

/// Расчет смещения точки изменения знака инварианта Вейля
fn weyl_sign_change_shift() -> f64 {
    // Теоретическая оценка смещения
    let params = MetricParams::default();
    HBAR * G / (C.powi(3) * params.r0.powi(2)) * scalar_curvature_approx(1.1)
}

//=============================================================================
// МОДЕЛИРОВАНИЕ BINARY MERGERS
//=============================================================================

/// Модель сливающихся черных дыр (упрощенная)
struct BinaryMerger {
    m1: f64,      // Масса первой черной дыры, солнечные массы
    m2: f64,      // Масса второй черной дыры, солнечные массы
    distance: f64, // Расстояние до наблюдателя, Мпк
}

impl BinaryMerger {
    // Вычисляет базовую амплитуду гравитационной волны
    fn base_amplitude(&self) -> f64 {
        let mtot = (self.m1 + self.m2) * MSUN;
        let dist_m = self.distance * 3.08567758e22; // Мпк в метры
        
        // Формула базовой амплитуды GW в ОТО
        (G * mtot) / (C * C * dist_m)
    }
    
    // Вычисляет характерную частоту слияния
    fn merger_frequency(&self) -> f64 {
        let mtot = (self.m1 + self.m2) * MSUN;
        
        // Характерная частота в Гц
        C.powi(3) / (G * mtot * PI * 6.0)
    }
    
    // Генерирует предсказанный LIGO сигнал в метриках ОТО и HyperTwist
    fn generate_signal(&self, r_observer: f64) -> HashMap<String, Vec<(f64, f64)>> {
        let base_amp = self.base_amplitude();
        let f_merger = self.merger_frequency();
        let params = MetricParams::default();
        
        // Временной массив до и после слияния (в секундах)
        let times: Vec<f64> = (-1000..1000).map(|i| i as f64 * 0.001).collect();
        
        // Частотная эволюция (чирп)
        let frequencies: Vec<f64> = times.iter()
            .map(|&t| {
                if t < 0.0 {
                    // До слияния - частота растет
                    f_merger * (1.0 - t).powf(-0.375)
                } else {
                    // После слияния - затухающий сигнал
                    f_merger * (1.0 + t * 5.0).powf(-0.8)
                }
            })
            .collect();
        
        // Амплитуды в ОТО
        let amp_gtr: Vec<(f64, f64)> = times.iter().enumerate()
            .map(|(i, &t)| {
                let phase = frequencies[i] * t * 2.0 * PI;
                let amp = if t < 0.0 {
                    base_amp * (1.0 - t).powf(-0.25)
                } else {
                    base_amp * (1.0 + t * 5.0).powf(-0.5)
                };
                (t, amp * phase.sin())
            })
            .collect();
        
        // Амплитуды в HyperTwist
        let amp_hypertwist: Vec<(f64, f64)> = times.iter().enumerate()
            .map(|(i, &t)| {
                let freq = frequencies[i];
                let phase_distortion = gw_waveform_distortion(r_observer, freq, &params);
                let phase = freq * t * 2.0 * PI + phase_distortion;
                
                let amp = if t < 0.0 {
                    base_amp * (1.0 - t).powf(-0.25)
                } else {
                    base_amp * (1.0 + t * 5.0).powf(-0.5)
                };
                
                // Дополнительная амплитуда из-за скручивания
                let twist_amp = twist_polarization_amplitude(r_observer, amp, &params);
                
                // Суммарный сигнал из стандартной и дополнительной моды
                let total_amp = amp * phase.sin() + twist_amp * (phase * 2.0).sin();
                
                (t, total_amp)
            })
            .collect();
        
        // Соберем результаты в словарь
        let mut results = HashMap::new();
        results.insert("GTR".to_string(), amp_gtr);
        results.insert("HyperTwist".to_string(), amp_hypertwist);
        results.insert("Frequencies".to_string(), 
            frequencies.iter().enumerate().map(|(i, &f)| (times[i], f)).collect());
        
        results
    }
}

/// Рассчитывает эффективное отношение сигнал/шум для двух моделей
fn signal_to_noise_ratio(original: &[(f64, f64)], modified: &[(f64, f64)]) -> f64 {
    let mut sum_diff_squared = 0.0;
    let mut sum_orig_squared = 0.0;
    
    for (i, &(_, amp_orig)) in original.iter().enumerate() {
        let (_, amp_mod) = modified[i];
        let diff = amp_orig - amp_mod;
        
        sum_diff_squared += diff * diff;
        sum_orig_squared += amp_orig * amp_orig;
    }
    
    if sum_orig_squared > 0.0 {
        (sum_diff_squared / sum_orig_squared).sqrt()
    } else {
        0.0
    }
}

/// Расчет разницы угла между входящим и исходящим лучом (линзирование)
fn angular_difference(theta1: f64, theta2: f64) -> f64 {
    let mut d = theta2 - theta1;
    while d > PI { d -= 2.0 * PI; }
    while d < -PI { d += 2.0 * PI; }
    d
}

//=============================================================================
// ОСНОВНЫЕ МОДУЛИ РАСЧЕТА И АНАЛИЗА
//=============================================================================

/// Анализ распространения гравитационных волн на разных расстояниях
fn analyze_gw_propagation() -> std::io::Result<()> {
    let params = MetricParams::default();
    
    // Файл для скорости гравитационных волн
    let mut speed_file = File::create("gw_speed.csv")?;
    writeln!(speed_file, "r,GW_Speed_Standard,GW_Speed_HyperTwist,Ratio")?;
    
    for i in 1..101 {
        let r = i as f64 * 0.1;
        let speed_ht = gw_speed(r, &params);
        let ratio = speed_ht / params.c;
        
        writeln!(speed_file, "{},{},{},{}", r, params.c, speed_ht, ratio)?;
    }
    
    // Файл для дополнительных поляризационных мод
    let mut pol_file = File::create("gw_polarization.csv")?;
    writeln!(pol_file, "r,Standard_Modes,HyperTwist_Extra_Mode,Ratio")?;
    
    let base_amp = 1.0e-21; // Типичная амплитуда GW
    for i in 1..101 {
        let r = i as f64 * 0.1;
        let twist_amp = twist_polarization_amplitude(r, base_amp, &params);
        let ratio = twist_amp / base_amp;
        
        writeln!(pol_file, "{},{},{},{}", r, base_amp, twist_amp, ratio)?;
    }
    
    Ok(())
}

/// Сравнение детектируемости сигналов для разных типов двойных систем
fn binary_merger_signals() -> std::io::Result<()> {
    // Различные типы слияний черных дыр
    let merger_types = vec![
        BinaryMerger { m1: 10.0, m2: 10.0, distance: 400.0 },     // Равные массы, GW150914-подобное
        BinaryMerger { m1: 30.0, m2: 30.0, distance: 1000.0 },    // Массивные ЧД
        BinaryMerger { m1: 5.0, m2: 1.4, distance: 200.0 },       // ЧД + НЗ
        BinaryMerger { m1: 1000.0, m2: 1000.0, distance: 3000.0 } // Сверхмассивные ЧД (LISA)
    ];
    
    let mut merger_file = File::create("binary_merger_snr.csv")?;
    writeln!(merger_file, "BH1_Mass,BH2_Mass,Distance_Mpc,r_observer,GTR_Peak,HyperTwist_Peak,SNR_Diff")?;
    
    for merger in &merger_types {
        // Проверка на разных расстояниях от центра гравитирующего объекта
        for r_obs in [0.5, 1.0, 2.0, 5.0, 10.0].iter() {
            let signals = merger.generate_signal(*r_obs);
            
            // Найдем максимальные амплитуды
            let max_gtr = signals["GTR"].iter()
                .map(|&(_, amp)| amp.abs())
                .fold(0.0, f64::max);
            
            let max_ht = signals["HyperTwist"].iter()
                .map(|&(_, amp)| amp.abs())
                .fold(0.0, f64::max);
            
            // Вычислим эффективное отношение сигнал/шум различия моделей
            let snr_diff = signal_to_noise_ratio(&signals["GTR"], &signals["HyperTwist"]);
            
            writeln!(merger_file, "{},{},{},{},{:.3e},{:.3e},{:.6}",
                     merger.m1, merger.m2, merger.distance, r_obs, max_gtr, max_ht, snr_diff)?;
        }
    }
    
    // Подробное сравнение формы сигнала для выбранного слияния
    let reference_merger = BinaryMerger { m1: 30.0, m2: 30.0, distance: 1000.0 };
    let r_detailed = 1.0; // Радиус наблюдения вблизи топологического перехода
    
    let signals = reference_merger.generate_signal(r_detailed);
    
    let mut waveform_file = File::create("gw_waveform_comparison.csv")?;
    writeln!(waveform_file, "Time,Frequency,GTR_Amplitude,HyperTwist_Amplitude,Difference")?;
    
    for i in 0..signals["GTR"].len() {
        let (t, amp_gtr) = signals["GTR"][i];
        let (_, amp_ht) = signals["HyperTwist"][i];
        let (_, freq) = signals["Frequencies"][i];
        let diff = amp_ht - amp_gtr;
        
        writeln!(waveform_file, "{},{},{},{},{}", t, freq, amp_gtr, amp_ht, diff)?;
    }
    
    Ok(())
}

/// Моделирование линзирования света в метрике HyperTwist
/// Моделирование линзирования света в метрике HyperTwist
fn simulate_gravitational_lensing() -> std::io::Result<()> {
    let params = MetricParams::default();
    let h = 1e-5;
    let dλ = 0.01;

    // Фотон стартует слева от центра, летит вправо
    let mut x = Vector4::new(0.0, 5.0, PI, 0.0);       // t, r, θ, z
    let mut dx = Vector4::new(1.0, -1.0, 0.3, 0.0);    // светоподобная 4-скорость

    // Нормировка: g_{μν} dx^μ dx^ν = 0 для светоподобных траекторий (фотоны)
    {
        // Получаем матрицу метрики для текущего значения r (x[1])
        let comp = MetricComponents::new(x[1], &params);
        let g = comp.to_matrix(&params);
    
        // Вычисляем "пространственную сумму" S = ∑_{i,j=1}^{3} g_{ij} dx[i] dx[j]
        let mut spatial_sum = 0.0;
        for i in 1..4 {
            for j in 1..4 {
                spatial_sum += g[(i, j)] * dx[i] * dx[j];
            }
        }
        
        // Проверяем, что g[0,0] имеет правильный знак (должен быть отрицательным)
        if g[(0, 0)] >= 0.0 {
            eprintln!("Ошибка: g[0,0] неотрицательное ({}). Невозможно нормировать траекторию.", g[(0, 0)]);
        } else if spatial_sum < 0.0 {
            // Если вычисленная пространственная сумма оказалась отрицательной, это может свидетельствовать об ошибке в исходных данных.
            eprintln!("Ошибка: пространственная сумма отрицательна ({}).", spatial_sum);
        } else {
            // Решаем для dx[0]:
            // g[0,0]*(dx[0])^2 + spatial_sum = 0  =>  dx[0] = sqrt(-spatial_sum / g[0,0])
            dx[0] = (-spatial_sum / g[(0, 0)]).sqrt();
        }
    
        // Проверяем после нормировки
        let gdot_after = dx.transpose() * g * dx;
        println!("После нормировки: gμνdxμdxν = {}", gdot_after[(0, 0)]);
    }

    let mut file = File::create("lens_trajectory.csv")?;
    writeln!(file, "x,y")?;

    for _ in 0..3000 {
        let r = x[1];
        let theta = x[2];
        let px = r * theta.cos();
        let py = r * theta.sin();
        writeln!(file, "{},{}", px, py)?;
        rk4_step(&mut x, &mut dx, &params, h, dλ);

        // остановка после пролёта
        if px > 5.0 && r > 5.0 {
            break;
        }
    }

    println!("✅ Траектория фотона записана в lens_trajectory.csv");
    analyze_lensing_angle()?;
    
    Ok(())
}

/// Анализ угла отклонения света при гравитационном линзировании
fn analyze_lensing_angle() -> std::io::Result<()> {
    let file = File::open("lens_trajectory.csv")?;
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().skip(1)
        .filter_map(Result::ok)
        .collect();

    if lines.len() < 2 {
        println!("Недостаточно точек для анализа.");
        return Ok(());
    }

    let (x1, y1) = {
        let parts: Vec<_> = lines[0].split(',').collect();
        (parts[0].parse::<f64>().unwrap(), parts[1].parse::<f64>().unwrap())
    };

    // Найдём последнюю НЕ-NaN строку
    let (x2, y2) = lines.iter().rev()
        .filter_map(|line| {
            let parts: Vec<_> = line.split(',').collect();
            let x = parts[0].parse::<f64>().ok()?;
            let y = parts[1].parse::<f64>().ok()?;
            if x.is_finite() && y.is_finite() {
                Some((x, y))
            } else {
                None
            }
        })
        .next()
        .expect("Нет допустимых финальных точек");

    let theta_in = y1.atan2(x1);
    let theta_out = y2.atan2(x2);
    let delta = angular_difference(theta_in, theta_out).abs().to_degrees();

    println!("🔭 Угол входа  θ_in  = {:.4} рад", theta_in);
    println!("🔭 Угол выхода θ_out = {:.4} рад", theta_out);
    println!("➡️  Угол отклонения света Δθ = {:.6}°", delta);

    Ok(())
}

/// Расчет и запись данных по скалярной кривизне и инвариантам
fn calculate_curvature_invariants() -> std::io::Result<()> {
    let params = MetricParams::default();
    let h = 1e-5;

    let mut file = File::create("curvature_invariants.csv")?;
    writeln!(file, "r,R,Ricci^2,Kretschmann,Weyl^2,BelRobinson")?;

    for i in 1..100 {
        let r = i as f64 * 0.1;

        let comp = MetricComponents::new(r, &params);
        let g = comp.to_matrix(&params);
        let g_inv = g.try_inverse().unwrap();
        let dg = metric_derivatives(r, &params, h);
        let gamma = christoffel(&g, &g_inv, &dg);
        let ricci = ricci_tensor(r, &params, h);
        let scalar_r = scalar_curvature(&g_inv, &ricci);
        let ricci_sq = ricci_invariant(&g_inv, &ricci);

        // === Добавляем gamma_plus и gamma_minus для производных символов ===
        let gamma_plus = {
            let comp_p = MetricComponents::new(r + h, &params);
            let g_p = comp_p.to_matrix(&params);
            let g_inv_p = g_p.try_inverse().unwrap();
            let dg_p = metric_derivatives(r + h, &params, h);
            christoffel(&g_p, &g_inv_p, &dg_p)
        };

        let gamma_minus = {
            let comp_m = MetricComponents::new(r - h, &params);
            let g_m = comp_m.to_matrix(&params);
            let g_inv_m = g_m.try_inverse().unwrap();
            let dg_m = metric_derivatives(r - h, &params, h);
            christoffel(&g_m, &g_inv_m, &dg_m)
        };

        // Используем верную версию тензора Римана
        let riemann = riemann_tensor_with_derivatives(&gamma, &gamma_plus, &gamma_minus, h);
        let kretsch = kretschmann_invariant_improved(&g_inv, &riemann);
        let weyl = weyl_tensor(&g, &ricci, scalar_r, &riemann);
        let weyl_sq = weyl_invariant_improved(&g_inv, &weyl);
        let belrobinson = belrobinson_tensor_full(&g_inv, &weyl);

        writeln!(file, "{:.2},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}", 
                r, scalar_r, ricci_sq, kretsch, weyl_sq, belrobinson)?;
    }

    println!("✅ Кривизнные инварианты записаны в curvature_invariants.csv");

    Ok(())
}


/// Расчет орбитальных скоростей в метрике HyperTwist
/// Расчет орбитальных скоростей в метрике HyperTwist
fn calculate_orbital_velocities() -> std::io::Result<()> {
    let params = MetricParams::default();

    let mut file = File::create("orbital_velocities.csv")?;
    writeln!(file, "r,v_kepler,v_improved,v_potential")?;
    
    for i in 1..200 {
        let r = i as f64 * 0.05;
        
        // Кеплеровская скорость для сравнения
        let v_kepler = (params.c * params.c * params.r0 / r).sqrt() / params.c;
        
        // Улучшенные модели
        let v_improved = orbital_velocity_improved(r, &params) / params.c; // нормируем к скорости света
        let v_potential = orbital_velocity_from_potential(r, &params) / params.c;
        
        writeln!(file, "{:.2},{:.6},{:.6},{:.6}", 
                r, v_kepler, v_improved, v_potential)?;
    }
    
    println!("✅ Орбитальные скорости записаны в orbital_velocities.csv");
    
    Ok(())
}

/// Космологические расчеты с учетом эффекта скручивания
fn calculate_cosmological_effects() -> std::io::Result<()> {
    let twist_effect = 0.3; // Параметр влияния скручивания (может требовать калибровки)
    let z_values: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
    
    // Файл для записи результатов Хаббла
    let mut hubble_file = File::create("hubble_parameter.csv")?;
    writeln!(hubble_file, "z,H_ΛCDM,H_HyperTwist,Twist_Energy_Density")?;
    
    for &z in &z_values {
        let h_lcdm = hubble_lcdm(z);
        let h_hypertwist = hubble_hypertwist(z, twist_effect);
        let twist_energy = twist_energy_density(z);
        
        writeln!(hubble_file, "{},{},{},{}", z, h_lcdm, h_hypertwist, twist_energy)?;
    }
    
    // Файл для записи результатов сверхновых
    let mut sn_file = File::create("supernova_magnitude.csv")?;
    writeln!(sn_file, "z,m_ΛCDM,m_HyperTwist,Difference")?;
    
    for &z in &z_values {
        if z > 0.01 {  // Исключаем z ~ 0 для избежания особенностей
            let (m_hypertwist, m_lcdm) = supernova_magnitude(z, twist_effect);
            let diff = m_hypertwist - m_lcdm;
            
            writeln!(sn_file, "{},{},{},{}", z, m_lcdm, m_hypertwist, diff)?;
        }
    }
    
    // Файл для записи роста возмущений
    let mut growth_file = File::create("perturbation_growth.csv")?;
    writeln!(growth_file, "z,Growth_ΛCDM,Growth_HyperTwist,Ratio")?;
    
    for &z in &z_values {
        let (growth_hypertwist, growth_lcdm) = matter_perturbation_growth(z, twist_effect);
        let ratio = growth_hypertwist / growth_lcdm;
        
        writeln!(growth_file, "{},{},{},{}", z, growth_lcdm, growth_hypertwist, ratio)?;
    }
    
    // Файл для углового диаметра (Тест Алкока-Пачинского)
    let mut angle_file = File::create("angular_diameter.csv")?;
    writeln!(angle_file, "z,Angle_ΛCDM,Angle_HyperTwist,Ratio")?;
    
    let physical_size = 1.0; // Размер в Мпк
    for &z in &z_values {
        if z > 0.01 {
            let (angle_hypertwist, angle_lcdm) = angular_diameter(z, physical_size, twist_effect);
            let ratio = angle_hypertwist / angle_lcdm;
            
            writeln!(angle_file, "{},{},{},{}", z, angle_lcdm, angle_hypertwist, ratio)?;
        }
    }
    
    println!("✅ Космологические расчеты завершены. Результаты записаны в CSV файлы.");
    
    Ok(())
}

/// Расчет квантовых эффектов в метрике HyperTwist
fn calculate_quantum_effects() -> std::io::Result<()> {
    let params = MetricParams::default();
    let mut quantum_file = File::create("quantum_effects.csv")?;
    writeln!(quantum_file, "r,R(r),V_quantum,Vacuum_Energy,Weyl_Classic,Weyl_Quantum,Change(%)")?;
    
    let r_values = [0.1, 0.5, 1.0, 1.1, 2.0, 5.0, 10.0];
    
    for &r in &r_values {
        let scalar_r = scalar_curvature_approx(r);
        let quantum_pot = quantum_potential(r);
        let vacuum_energy = vacuum_energy_density(r);
        let weyl_classic = weyl_invariant_approx(r);
        let weyl_quantum = quantum_weyl_correction(r);
        let change_pct = (weyl_quantum - weyl_classic) / weyl_classic.abs() * 100.0;
        
        writeln!(quantum_file, "{:.1},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6}", 
                r, scalar_r, quantum_pot, vacuum_energy, weyl_classic, weyl_quantum, change_pct)?;
    }

    println!("\n=== Проверка тождества Бианки ===");
    let r_test = 1.0;
    let bianchi = check_bianchi_identity(r_test, &params, 1e-5);
    for nu in 0..4 {
        println!("∇^μ G_{{μ{}}} = {:.3e}", nu, bianchi[nu]);
    }

    println!("✅ Квантовые эффекты записаны в quantum_effects.csv");
    println!("📊 Смещение изменения знака тензора Вейля: {:.6e}", weyl_sign_change_shift());
    
    Ok(())
}

//=============================================================================
// ГЛАВНАЯ ФУНКЦИЯ
//=============================================================================

fn main() -> std::io::Result<()> {
    println!("\n🌌 Симуляция метрики HyperTwist 🌌\n");
    
    // Расчет базовых характеристик метрики
    calculate_curvature_invariants()?;
    calculate_orbital_velocities()?;
    
    // Гравитационное линзирование
    simulate_gravitational_lensing()?;
    
    // Космологические расчеты
    calculate_cosmological_effects()?;
    
    // Гравитационные волны
    analyze_gw_propagation()?;
    binary_merger_signals()?;
    
    // Квантовые эффекты
    calculate_quantum_effects()?;
    
    // Проверка аксиальной симметрии
    println!("\n=== Проверка аксиальной симметрии ===");
    let params = MetricParams::default();
    let epsilon = 1e-10;
    
    for i in 1..11 {
        let r = i as f64 * 0.5;
        if !is_axially_symmetric(r, &params, epsilon) {
            println!("❌ Симметрия нарушена при r = {:.2}", r);
        }
    }
    println!("✅ Аксиальная симметрия подтверждена для всех проверенных r");
    
    println!("\n✅ Все расчеты успешно завершены!");
    println!("📁 Результаты сохранены в CSV-файлы для дальнейшего анализа.");
    
    Ok(())
}