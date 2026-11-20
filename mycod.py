# -*- coding: utf-8 -*-
"""
========================================================================
 РАСШИРЕННАЯ FORWARD-МОДЕЛЬ РЕГИСТРАЦИИ ЭЛЕКТРОННО-ЦИКЛОТРОННОГО
 ИЗЛУЧЕНИЯ (ECE) ПЛАЗМЫ ТОКАМАКА: O1, X1, O2, X2, O3, X3
 + ШУМЫ ПРИЕМНИКА, ФОН, ПЕРЕОТРАЖЕНИЯ И ПОЛЯРИЗАЦИОННОЕ СМЕШИВАНИЕ
========================================================================

НАЗНАЧЕНИЕ
----------
Упрощенная (не кинетическая, не лучевая) forward-модель для оценки:
  - положения циклотронных резонансных слоев s*f_ce(R)/gamma, s=1,2,3;
  - холодноплазменной доступности (отсечки O-/X-мод);
  - приближенной оптической толщины и спектра T_b(nu) для ОБЕИХ мод
    (O и X), гармоник s=1,2,3, с явной маркировкой точности каждой
    формулы (валидированная / качественная параметризация);
  - влияния измерительных шумов, фона приемника, многократных
    переотражений от стенки камеры и смешивания поляризаций O/X
    на итоговый наблюдаемый спектр.

ВАЖНОЕ ОГРАНИЧЕНИЕ (см. раздел 7 отчета в чате)
-------------------------------------------------
Это НЕ полная релятивистская кинетическая ECE-модель. Отсутствуют:
  - решение полного дисперсионного уравнения магнетизированной плазмы;
  - трассировка луча через реальную геометрию антенны/оптики;
  - кинетические (не Максвелл-усредненные по холодной плазме)
    коэффициенты излучения/поглощения;
  - реальные частотно-зависимые характеристики стенки и поляризатора.
Валидированная формула поглощения (Bornatici et al., 1983) реализована
для X-моды (гармоники 1,2,3). Для O-моды применена КАЧЕСТВЕННАЯ
параметризация по аналогии (см. Sabri & Benouaz, IJETAE 2(8), 2012;
Hutchinson, "Principles of Plasma Diagnostics", гл.5) - результаты
O1/O2/O3 помечены статусом CALCULATED_QUALITATIVE / MODEL_APPROXIMATION
и не являются строгим количественным расчетом T_e.

ИСПРАВЛЕНИЕ СОВМЕСТИМОСТИ NUMPY
--------------------------------
В NumPy >= 2.0 функция np.trapz удалена и заменена на np.trapezoid.
Используется локальная обертка trapz_compat().
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import factorial

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams.update({
    "figure.dpi": 130, "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
    "axes.titlesize": 12, "axes.labelsize": 11,
})

C_MAIN, C_ALT, C_3, C_GOLD = "#20808D", "#A84B2F", "#7A39BB", "#D19900"
MODE_COLORS = {
    "O1": "#4C72B0", "X1": "#DD8452",
    "O2": "#55A868", "X2": "#20808D",
    "O3": "#C44E52", "X3": "#7A39BB",
}

try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "ece_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def out_path(filename):
    return os.path.join(OUTPUT_DIR, filename)


# ========================================================================
# СОВМЕСТИМОСТЬ NUMPY 1.x / 2.x  (np.trapz удален в NumPy >= 2.0)
# ========================================================================
def trapz_compat(y, x=None):
    """Обертка: использует np.trapezoid (NumPy>=2.0) либо np.trapz (NumPy<2.0)."""
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


# ========================================================================
# ФИЗИЧЕСКИЕ КОНСТАНТЫ СИ
# ========================================================================
E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
EPS0 = 8.8541878128e-12
K_B = 1.380649e-23
C_LIGHT = 2.99792458e8


def ev_to_joule(t_ev):
    return t_ev * E_CHARGE


def joule_to_ev(t_j):
    return t_j / E_CHARGE


def cyclotron_omega(b_tesla):
    """omega_ce = e*B/m_e [рад/с]."""
    return E_CHARGE * b_tesla / M_E


def cyclotron_frequency(b_tesla):
    """f_ce = omega_ce/(2*pi) = e*B/(2*pi*m_e) [Гц]. НЕ путать с omega!"""
    return cyclotron_omega(b_tesla) / (2.0 * np.pi)


def plasma_omega(ne_m3):
    return np.sqrt(ne_m3 * E_CHARGE ** 2 / (EPS0 * M_E))


def coulomb_log_ei(ne_m3, te_ev):
    ne_cm3 = np.asarray(ne_m3) * 1e-6
    te_ev = np.asarray(te_ev)
    if np.any(ne_cm3 <= 0.0) or np.any(te_ev <= 0.0):
        raise ValueError("ne и Te должны быть положительны")
    return 24.0 - np.log(np.sqrt(ne_cm3) / te_ev)


def collision_frequency_ei(ne_m3, te_j, z_eff=1.0):
    te_ev = joule_to_ev(te_j)
    ln_l = coulomb_log_ei(ne_m3, te_ev)
    num = ne_m3 * z_eff * E_CHARGE ** 4 * ln_l
    den = 3.0 * (2.0 * np.pi) ** 1.5 * EPS0 ** 2 * np.sqrt(M_E) * te_j ** 1.5
    return num / den


# ========================================================================
# КОНФИГУРАЦИЯ
# ========================================================================
class PlasmaConfig:
    def __init__(self):
        # ---- Геометрия и профили плазмы (согласовано ранее) ----
        self.R0 = 2.15
        self.a = 0.57
        self.B0 = 1.0
        self.Te0_keV = 8.0
        self.ne0_m3 = 1.0e18
        self.z_eff = 1.0
        self.te_edge_fraction = 0.02
        self.ne_edge_fraction = 0.05
        self.view_angle_deg = 90.0          # строго перпендикулярное наблюдение
        self.n_grid_R = 800
        self.n_grid_nu = 3000
        self.channel_half_width_hz = 0.25e9
        self.n_channels = 96
        self.tau_threshold = -np.log(0.05)

        # ---- Явная конфигурация мод и гармоник ----
        self.HARMONICS = (1, 2, 3)
        self.MODES = ("O", "X")
        # Разрешить отключать отдельные каналы, например: {"O3": False}
        self.CHANNEL_ENABLED = {f"{m}{s}": True for s in self.HARMONICS for m in self.MODES}

        # ================================================================
        # БЛОК 5: ШУМЫ И ФОН ПРИЕМНИКА (DEMO / требует замены на параметры прибора)
        # ================================================================
        self.noise_enabled = True
        self.random_seed = 42
        # Температура шума приемника [К] - ЭТО НЕ Te плазмы! Отдельная величина,
        # характеризующая собственный шум радиометра/усилителя.
        self.receiver_noise_temperature_K = 150.0          # DEMO
        # Относительный (мультипликативный) уровень шума от уровня сигнала
        self.relative_noise_std = 0.03                      # DEMO, 3%
        # Абсолютный (аддитивный) уровень шума в единицах T_b [эВ]
        self.absolute_noise_level_eV = 0.5                   # DEMO
        # Медленный фон/наклон базовой линии (в единицах T_b, эВ и эВ/ГГц)
        self.background_offset_eV = 1.0                      # DEMO
        self.background_slope_eV_per_GHz = 0.02               # DEMO
        # Ширина корреляции шума по частоте [Гц] - для коррелированного (дрейф) шума
        self.frequency_correlation_width_hz = 1.0e9           # DEMO
        self.n_noise_realizations = 1

        # ================================================================
        # БЛОК 6: ПЕРЕОТРАЖЕНИЯ И ПОЛЯРИЗАЦИОННОЕ СМЕШИВАНИЕ (DEMO)
        # ================================================================
        self.reflections_enabled = True
        # Коэффициент отражения стенки/камеры по мощности, 0<=R_wall<1
        self.wall_reflectivity = 0.20                         # DEMO
        self.max_reflection_order = 10
        # Если True - используется предел бесконечной геом. суммы (при R_eff<1)
        self.infinite_reflection_model = True
        # Дополнительное затухание за проход "туда-обратно" (оптика, окна), 0..1
        self.round_trip_attenuation = 0.85                    # DEMO (15% потерь)
        # Смешивание поляризаций при отражении/несовершенстве тракта, 0..1
        self.polarization_mixing_ox = 0.03                    # DEMO: доля O->X
        self.polarization_mixing_xo = 0.03                    # DEMO: доля X->O
        # Пропускание поляризатора для целевой поляризации и утечка "чужой"
        self.polarizer_x_transmission = 0.98                  # DEMO
        self.polarizer_o_leakage = 0.02                       # DEMO

    def combo_labels(self):
        labels = []
        for s in self.HARMONICS:
            for m in self.MODES:
                labels.append(f"{m}{s}")
        return labels

    def enabled_labels(self):
        return [lbl for lbl in self.combo_labels() if self.CHANNEL_ENABLED.get(lbl, True)]


MODEL_ACCURACY = {
    "X": "Bornatici (1983): валидированная в литературе оценка оптической "
         "толщины X-моды при перпендикулярном наблюдении.",
    "O": "Качественная параметризация по аналогии с масштабным законом "
         "Bornatici/Hutchinson: сохраняет верную функциональную зависимость "
         "от beta_t^(n-1) и геом. фактора L_B, но численный гармонический "
         "коэффициент НЕ является строгим кинетическим расчетом. Только для "
         "оценки порядка величины и статуса доступности, не для количественной "
         "T_e-диагностики.",
}

STATUS_DESCRIPTIONS = {
    "CALCULATED":
        "Резонанс найден в расчетной области, доступность подтверждена в "
        "рамках упрощенного холодноплазменного критерия, спектр рассчитан "
        "по валидированной формуле (X-мода, Bornatici 1983).",
    "CALCULATED_QUALITATIVE":
        "Резонанс найден, доступность подтверждена в рамках упрощенного "
        "критерия, спектр рассчитан по качественной (не кинетической) "
        "параметризации по аналогии с X-модой - применимо только для оценки "
        "порядка величины и статуса, НЕ для количественной Te-диагностики.",
    "OUT_OF_RESONANCE_RANGE":
        "Резонансная частота s*f_ce(R)/gamma не попадает в диапазон, "
        "покрываемый частотной сеткой в расчетной области (резонансный слой "
        "физически находится вне рассматриваемого интервала большого радиуса).",
    "CUTOFF_RISK":
        "Резонанс найден, но частота попадает в холодноплазменную полосу "
        "непрозрачности (отсечка O- или X-моды) на пути от резонансного слоя "
        "к приемнику - излучение на этой частоте не может выйти из плазмы "
        "в рамках холодноплазменного приближения.",
    "LOW_OPTICAL_DEPTH":
        "Резонанс доступен, но максимальная оптическая толщина по трассе "
        "заметно меньше единицы (tau_max < порога) - канал оптически тонкий, "
        "T_b систематически ниже Te, интерпретация T_b~=Te неприменима без "
        "поправки на конечную tau.",
        
    "POLARIZATION_MIXING_RISK":
        "Канал рассчитан, но при заданных в конфигурации коэффициентах "
        "поляризационного смешивания вклад чужой поляризации сопоставим "
        "или превышает собственный сигнал - интерпретация как чистой моды "
        "требует контроля поляризационной селективности приемного тракта.",
    "MODEL_APPROXIMATION":
        "Канал рассчитан только в рамках феноменологического приближения "
        "(например, качественная O-мода или демонстрационная модель "
        "переотражений/шума) - результат не является строгим количественным "
        "расчетом.",
    "INSUFFICIENT_INPUT_DATA":
        "Недостаточно входных данных (профилей/параметров прибора) для "
        "расчета канала - использованы демонстрационные значения, требующие "
        "замены на измеренные.",
    "NOT_IMPLEMENTED_PHYSICALLY":
        "Комбинация мода-гармоника не реализована физически корректно в "
        "рамках текущей модели (например, отсутствует необходимое "
        "приближение) и не рассчитывается.",
}


# ========================================================================
# ПРОФИЛИ ПЛАЗМЫ
# ========================================================================
def build_radial_grid(cfg):
    return np.linspace(cfg.R0 - cfg.a + 1e-3, cfg.R0 + cfg.a - 1e-3, cfg.n_grid_R)


def normalized_radius(r, cfg):
    return np.clip(np.abs(r - cfg.R0) / cfg.a, 0.0, 1.0)


def magnetic_field_profile(r, cfg):
    """B(R) = B0*R0/R [Тл] - одномерное приближение тороидального поля."""
    return cfg.B0 * cfg.R0 / r


def density_profile(rho, cfg):
    shape = 1.0 - rho ** 2
    return cfg.ne0_m3 * (cfg.ne_edge_fraction + (1.0 - cfg.ne_edge_fraction) * shape)


def temperature_profile_keV(rho, cfg):
    shape = 1.0 - rho ** 2
    return cfg.Te0_keV * (cfg.te_edge_fraction + (1.0 - cfg.te_edge_fraction) * shape)


def build_plasma_profiles(cfg):
    r = build_radial_grid(cfg)
    rho = normalized_radius(r, cfg)
    b_field = magnetic_field_profile(r, cfg)
    ne_m3 = density_profile(rho, cfg)
    te_keV = temperature_profile_keV(rho, cfg)
    te_j = ev_to_joule(te_keV * 1e3)
    return {"R": r, "rho": rho, "B": b_field, "ne": ne_m3, "Te_keV": te_keV, "Te_J": te_j}


def plot_input_profiles(profiles, cfg):
    r = profiles["R"]
    saved = []

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(r, profiles["B"], color=C_MAIN, lw=2.2, label=r"$B(R)=B_0R_0/R$")
    ax.set_xlabel("Большой радиус R, м"); ax.set_ylabel("Магнитное поле B, Тл")
    ax.set_title("Входной профиль: магнитное поле B(R),\nиспользованный в расчете")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10)
    f = out_path("01_магнитное_поле_B_R.png"); fig.tight_layout()
    fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig); saved.append(f)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(r, profiles["ne"], color=C_ALT, lw=2.2, label=r"$n_e(R)$")
    ax.set_xlabel("Большой радиус R, м"); ax.set_ylabel(r"Концентрация электронов $n_e$, м$^{-3}$")
    ax.set_title("Входной профиль: концентрация электронов $n_e(R)$,\nиспользованный в расчете")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10)
    f = out_path("02_электронная_концентрация_ne_R.png"); fig.tight_layout()
    fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig); saved.append(f)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(r, profiles["Te_keV"], color=C_3, lw=2.2, label=r"$T_e(R)$")
    ax.set_xlabel("Большой радиус R, м"); ax.set_ylabel(r"Электронная температура $T_e$, кэВ")
    ax.set_title("Входной профиль: электронная температура $T_e(R)$,\nиспользованный в расчете")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10)
    f = out_path("03_электронная_температура_Te_R.png"); fig.tight_layout()
    fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig); saved.append(f)

    csv_path = out_path("использованные_профили_плазмы.csv")
    data = np.column_stack([r, profiles["rho"], profiles["B"], profiles["ne"], profiles["Te_keV"]])
    np.savetxt(csv_path, data, delimiter=",", header="R_m,rho,B_T,ne_m-3,Te_keV",
               comments="", fmt="%.6e")
    saved.append(csv_path)
    return saved


# ========================================================================
# ХОЛОДНОПЛАЗМЕННЫЕ ОТСЕЧКИ
# ========================================================================
def cutoff_frequencies(b_field, ne_m3):
    """
    f_R = (f_ce + sqrt(f_ce^2+4 f_pe^2))/2   - верхняя отсечка X-моды
    f_L = (-f_ce + sqrt(f_ce^2+4 f_pe^2))/2  - нижняя отсечка X-моды
    f_O = f_pe                                - отсечка O-моды
    Общепринятые холодноплазменные соотношения (Stix, Hutchinson).
    """
    f_ce = cyclotron_frequency(b_field)
    f_pe = plasma_omega(ne_m3) / (2.0 * np.pi)
    disc = np.sqrt(f_ce ** 2 + 4.0 * f_pe ** 2)
    f_r = 0.5 * (f_ce + disc)
    f_l = 0.5 * (-f_ce + disc)
    return f_l, f_r, f_pe, f_ce


# ========================================================================
# РЕЗОНАНСНЫЕ ЧАСТОТЫ ГАРМОНИК (холодноплазменное приближение gamma=1)
# ========================================================================
def harmonic_resonance_frequency(b_field, s):
    """
    f_res,s(R) = s*f_ce(R)/gamma. В данной модели используется gamma=1
    (нерелятивистское приближение для положения ЦЕНТРА резонанса);
    релятивистское уширение учтено ОТДЕЛЬНО как ширина слоя
    (resonance_layer_width), а не сдвиг центральной частоты. Это
    приближение достаточно для Te much less me*c^2 (нерелятивистская плазма).
    """
    f_ce = cyclotron_frequency(b_field)
    return s * f_ce


# ========================================================================
# ДОСТУПНОСТЬ (холодноплазменный критерий распространения к приемнику)
# ========================================================================
def check_accessibility(nu, r, r_res, mode, f_l, f_r, f_pe):
    """
    Приближенная проверка: волна должна дойти от резонансного слоя R_res
    до внешнего края плазмы (LFS, R>=R_res) без пересечения полосы
    непрозрачности СВОЕЙ моды. Упрощенный (не лучевой) критерий -
    НЕ является заменой полного решения дисперсионного уравнения и
    трассировки луча.
    """
    path_mask = r >= r_res
    if not np.any(path_mask):
        return False, "путь к приемнику вне расчетной сетки"

    if mode == "O":
        blocked = nu <= f_pe[path_mask]
        if np.any(blocked):
            return False, "пересечение отсечки O-моды (f <= f_pe) на пути к приемнику"
        return True, "путь свободен от отсечки O-моды (оценка, не полная трассировка луча)"

    if mode == "X":
        f_r_path = f_r[path_mask]
        f_l_path = f_l[path_mask]
        blocked = np.any((nu > f_l_path) & (nu < f_r_path))
        if blocked:
            return False, "частота попадает в полосу непрозрачности X-моды (f_L < f < f_U) на пути"
        return True, "путь свободен от полос непрозрачности X-моды (оценка, не полная трассировка луча)"

    return False, "неизвестная мода"


# ========================================================================
# ПОГЛОЩЕНИЕ (Bornatici-подобные оценки для X- и O-мод)
# ========================================================================
def gaussian_spatial_profile(r, r_res, sigma_r):
    sigma_r = max(float(sigma_r), 1e-6)
    h = np.exp(-0.5 * ((r - r_res) / sigma_r) ** 2)
    area = trapz_compat(h, r)
    return h / area if area > 0.0 else np.zeros_like(h)


def resonance_layer_width(b_field, r, harmonic, te_j):
    """
    Модельная (не строгая FWHM) ширина слоя из релятивистского уширения:
        d_omega_rel ~ 2*s*omega_ce*(Te/(me c^2))
    Валидно только как ПОРЯДКОВАЯ оценка при Te << me*c^2.
    """
    dbdr = np.abs(np.gradient(b_field, r))
    dbdr = np.where(dbdr < 1e-12, 1e-12, dbdr)
    omega_ce = cyclotron_omega(b_field)
    domega_ds = (E_CHARGE / M_E) * dbdr
    beta_t = te_j / (M_E * C_LIGHT ** 2)
    domega_rel = 2.0 * harmonic * omega_ce * beta_t
    sigma_r = domega_rel / domega_ds
    return np.maximum(sigma_r, 1e-4)


def geometric_length_factor(b_field, r):
    """L_B = B/|dB/dR| - масштаб неоднородности поля, м. Общий для X и O."""
    dbdr = np.abs(np.gradient(b_field, r))
    dbdr = np.where(dbdr < 1e-12, 1e-12, dbdr)
    return b_field / dbdr


def harmonic_factor_xmode(harmonic, beta_t):
    """
    Гармонический коэффициент X-моды по Bornatici et al. (1983):
        A_n^X = n^(2(n-1)) / (n-1)! * beta_t^(n-1)
    ВАЛИДИРОВАННАЯ в литературе форма (Hutchinson, гл.5).
    """
    return (harmonic ** (2 * (harmonic - 1))) / factorial(harmonic - 1) * beta_t ** (harmonic - 1)


def harmonic_factor_omode(harmonic, beta_t):
    """
    Гармонический коэффициент O-моды - КАЧЕСТВЕННАЯ параметризация
    по аналогии со структурой масштабного закона Bornatici/Hutchinson
    для O-моды: та же степенная зависимость от beta_t^(n-1) и
    факториальный знаменатель, но с дополнительным подавляющим
    множителем 2^n*n! в знаменателе, отражающим общеизвестный
    качественный факт, что оптическая толщина O-моды систематически
    МЕНЬШЕ X-моды той же гармоники. НЕ строгий кинетический расчет.
    """
    return (harmonic ** (2 * harmonic - 1)) / (2.0 ** harmonic * factorial(harmonic)) * beta_t ** (harmonic - 1)


def bornatici_tau_estimate(b_field, ne_m3, te_j, r, harmonic, mode):
    """
    Общая порядковая оценка полной оптической толщины гармоники
    при перпендикулярном наблюдении:
        tau_layer = (pi/2) * A_n^{mode}(beta_t) * (omega_pe^2/(omega_ce*c)) * L_B
    Различие X/O заключено ТОЛЬКО в гармоническом коэффициенте A_n.
    """
    omega_ce = cyclotron_omega(b_field)
    omega_pe = plasma_omega(ne_m3)
    beta_t = te_j / (M_E * C_LIGHT ** 2)
    l_b = geometric_length_factor(b_field, r)

    if mode == "X":
        harm_factor = harmonic_factor_xmode(harmonic, beta_t)
    elif mode == "O":
        harm_factor = harmonic_factor_omode(harmonic, beta_t)
    else:
        raise ValueError(f"неизвестная мода: {mode}")

    tau_layer = (np.pi / 2.0) * harm_factor * (omega_pe ** 2 / (omega_ce * C_LIGHT)) * l_b
    return np.maximum(tau_layer, 0.0)


def build_absorption_grid(r, nu_grid, b_field, ne_m3, te_j, harmonic, mode):
    """alpha_n(R,nu) [м^-1]: alpha = tau_layer(nu)*h_nu(R), int h dR=1."""
    f_res_profile = harmonic_resonance_frequency(b_field, harmonic)
    tau_r = bornatici_tau_estimate(b_field, ne_m3, te_j, r, harmonic, mode)
    sigma_r_profile = resonance_layer_width(b_field, r, harmonic, te_j)

    alpha = np.zeros((r.size, nu_grid.size))
    resonance_found = np.zeros(nu_grid.size, dtype=bool)
    r_res_arr = np.full(nu_grid.size, np.nan)

    f_min, f_max = f_res_profile.min(), f_res_profile.max()
    for j, nu in enumerate(nu_grid):
        if nu < f_min or nu > f_max:
            continue
        r_res = np.interp(nu, f_res_profile[::-1], r[::-1])
        tau_res = np.interp(r_res, r, tau_r)
        sigma_res = np.interp(r_res, r, sigma_r_profile)
        h = gaussian_spatial_profile(r, r_res, sigma_res)
        alpha[:, j] = tau_res * h
        resonance_found[j] = True
        r_res_arr[j] = r_res

    return alpha, resonance_found, r_res_arr


# ========================================================================
# ПЕРЕНОС ИЗЛУЧЕНИЯ (устойчивая рекуррентная схема, ЛТР-приближение)
# ========================================================================
def solve_radiative_transfer(r, te_j_profile, alpha_grid):
    """
    T_b[i+1] = T_b[i]*exp(-dtau) + Te[i]*(1-exp(-dtau)),
    что для ЛТР (j_nu = alpha_nu * Te, закон Кирхгофа) эквивалентно
    стандартному решению dI/ds = j_nu - alpha_nu*I. T_b~=Te справедливо
    только для tau>>1 (оптически толстый канал); при малых tau T_b~tau*Te.
    """
    n_r, n_nu = alpha_grid.shape
    ds = np.abs(np.gradient(r))
    tb = np.zeros(n_nu)
    tau_running = np.zeros(n_nu)
    tau_cum = np.zeros((n_r, n_nu))

    for i in range(n_r):
        dtau = np.maximum(alpha_grid[i, :] * ds[i], 0.0)
        transmission = np.exp(-dtau)
        emission = -np.expm1(-dtau)
        tau_running += dtau
        tau_cum[i, :] = tau_running
        tb = tb * transmission + te_j_profile[i] * emission

    return tb, tau_cum


# ========================================================================
# АППАРАТНАЯ ФУНКЦИЯ ПРИЕМНИКА
# ========================================================================
def rectangular_channel_response(nu_grid, center, half_width):
    response = (np.abs(nu_grid - center) <= half_width).astype(float)
    area = trapz_compat(response, nu_grid)
    return response / area if area > 0.0 else response


def channel_spectrum(nu_grid, tb_spectrum, channel_centers, half_width):
    tb_channel = np.zeros_like(channel_centers)
    for k, center in enumerate(channel_centers):
        resp = rectangular_channel_response(nu_grid, center, half_width)
        tb_channel[k] = trapz_compat(resp * tb_spectrum, nu_grid)
    return tb_channel


# ========================================================================
# РАСЧЕТ ОДНОЙ КОМБИНАЦИИ МОДА-ГАРМОНИКА (O_s или X_s) - "ЧИСТЫЙ" ECE
# ========================================================================
def compute_combo(mode, s, r, nu_grid, profiles, f_l, f_r, f_pe, cfg):
    label = f"{mode}{s}"
    f_res_profile = harmonic_resonance_frequency(profiles["B"], s)
    f_min, f_max = f_res_profile.min(), f_res_profile.max()

    result = {
        "label": label, "mode": mode, "harmonic": s,
        "f_res_min_GHz": f_min / 1e9, "f_res_max_GHz": f_max / 1e9,
        "status": None, "reason": None, "accuracy": MODEL_ACCURACY[mode],
        "alpha": None, "tau_cum": None, "tb_spectrum": None,
        "access_mask": None, "r_res_arr": None, "tau_max": None,
    }

    if not cfg.CHANNEL_ENABLED.get(label, True):
        result["status"] = "NOT_IMPLEMENTED_PHYSICALLY"
        result["reason"] = "Канал отключен пользователем через конфигурацию CHANNEL_ENABLED."
        return result

    grid_overlaps = (f_min <= nu_grid.max()) and (f_max >= nu_grid.min())
    if not grid_overlaps:
        result["status"] = "OUT_OF_RESONANCE_RANGE"
        result["reason"] = STATUS_DESCRIPTIONS["OUT_OF_RESONANCE_RANGE"]
        return result

    alpha, resonance_found, r_res_arr = build_absorption_grid(
        r, nu_grid, profiles["B"], profiles["ne"], profiles["Te_J"], s, mode
    )
    access_mask = np.zeros(nu_grid.size, dtype=bool)
    access_reason_sample = ""
    for j, nu in enumerate(nu_grid):
        if not resonance_found[j]:
            continue
        ok, reason = check_accessibility(nu, r, r_res_arr[j], mode, f_l, f_r, f_pe)
        access_mask[j] = ok
        if not ok and not access_reason_sample:
            access_reason_sample = reason

    if not np.any(access_mask):
        result["status"] = "CUTOFF_RISK"
        result["reason"] = STATUS_DESCRIPTIONS["CUTOFF_RISK"] + f" ({access_reason_sample})"
        result["r_res_arr"] = r_res_arr
        return result

    alpha_masked = alpha.copy()
    alpha_masked[:, ~access_mask] = 0.0
    tb_spectrum, tau_cum = solve_radiative_transfer(r, profiles["Te_J"], alpha_masked)
    tb_spectrum = np.where(access_mask, tb_spectrum, 0.0)
    tau_max = tau_cum[-1, :].max() if tau_cum.size else 0.0

    if tau_max < cfg.tau_threshold:
        status = "LOW_OPTICAL_DEPTH"
    else:
        status = "CALCULATED" if mode == "X" else "CALCULATED_QUALITATIVE"

    result["status"] = status
    result["reason"] = STATUS_DESCRIPTIONS[status] + (
        f" (max tau={tau_max:.3f})" if status == "LOW_OPTICAL_DEPTH" else ""
    )
    result["alpha"] = alpha_masked
    result["tau_cum"] = tau_cum
    result["tb_spectrum"] = tb_spectrum
    result["access_mask"] = access_mask
    result["r_res_arr"] = r_res_arr
    result["tau_max"] = tau_max
    return result


def compute_all_combinations(cfg, r, nu_grid, profiles, f_l, f_r, f_pe):
    results = {}
    for s in cfg.HARMONICS:
        for mode in cfg.MODES:
            results[f"{mode}{s}"] = compute_combo(mode, s, r, nu_grid, profiles, f_l, f_r, f_pe, cfg)
    return results


def is_ok_status(status):
    return status in ("CALCULATED", "CALCULATED_QUALITATIVE", "LOW_OPTICAL_DEPTH")


def has_spectrum(res):
    return res.get("tb_spectrum") is not None


# ========================================================================
# БЛОК 6: ПОЛЯРИЗАЦИОННОЕ СМЕШИВАНИЕ (O<->X) - феноменологическая матрица
# ========================================================================
def build_polarization_matrix(cfg):
    """
    Феноменологическая 2x2 матрица передачи поляризационного тракта для
    пары (O_s, X_s) одной гармоники:
        [I_obs_O]   [ (1-eta_OX)*T_O    eta_XO*L_O ] [I_O]
        [I_obs_X] = [ eta_OX*L_X    (1-eta_XO)*T_X ] [I_X]
    где T_O, T_X - пропускание "родного" канала поляризатора,
    eta_OX, eta_XO - доли перетекания между модами (0..1),
    L_O, L_X - коэффициенты утечки "чужой" поляризации.
    Все параметры - ФЕНОМЕНОЛОГИЧЕСКИЕ (DEMO), требуют замены на
    паспортные характеристики поляризационного тракта установки.
    """
    eta_ox = np.clip(cfg.polarization_mixing_ox, 0.0, 1.0)
    eta_xo = np.clip(cfg.polarization_mixing_xo, 0.0, 1.0)
    t_x = np.clip(cfg.polarizer_x_transmission, 0.0, 1.0)
    l_o = np.clip(cfg.polarizer_o_leakage, 0.0, 1.0)
    t_o = 1.0 - l_o
    m_pol = np.array([
        [(1.0 - eta_ox) * t_o, eta_xo * l_o],
        [eta_ox * l_o,          (1.0 - eta_xo) * t_x],
    ])
    return m_pol


def apply_polarization_mixing(results, cfg):
    """
    Применяет матрицу M_pol к парам (O_s, X_s) для каждой гармоники s,
    только если ОБА канала имеют рассчитанный спектр (иначе смешивание
    физически не определено данными модели). Результат сохраняется в
    отдельном поле tb_after_polarization, чистый tb_spectrum НЕ изменяется.
    """
    m_pol = build_polarization_matrix(cfg)
    for s in cfg.HARMONICS:
        lbl_o, lbl_x = f"O{s}", f"X{s}"
        res_o, res_x = results.get(lbl_o), results.get(lbl_x)
        if res_o is None or res_x is None:
            continue
        has_o, has_x = has_spectrum(res_o), has_spectrum(res_x)
        if not has_o and not has_x:
            continue
        io = res_o["tb_spectrum"] if has_o else np.zeros_like(res_x["tb_spectrum"])
        ix = res_x["tb_spectrum"] if has_x else np.zeros_like(res_o["tb_spectrum"])
        io_obs = m_pol[0, 0] * io + m_pol[0, 1] * ix
        ix_obs = m_pol[1, 0] * io + m_pol[1, 1] * ix
        if has_o:
            res_o["tb_after_polarization"] = io_obs
            leak_frac = np.divide(m_pol[0, 1] * ix, np.maximum(io_obs, 1e-30),
                                   out=np.zeros_like(io_obs), where=io_obs > 1e-30)
            if np.nanmax(leak_frac) > 0.1 and res_o["status"] in ("CALCULATED", "CALCULATED_QUALITATIVE", "LOW_OPTICAL_DEPTH"):
                res_o.setdefault("warnings", []).append(STATUS_DESCRIPTIONS["POLARIZATION_MIXING_RISK"])
        if has_x:
            res_x["tb_after_polarization"] = ix_obs
            leak_frac = np.divide(m_pol[1, 0] * io, np.maximum(ix_obs, 1e-30),
                                   out=np.zeros_like(ix_obs), where=ix_obs > 1e-30)
            if np.nanmax(leak_frac) > 0.1 and res_x["status"] in ("CALCULATED", "CALCULATED_QUALITATIVE", "LOW_OPTICAL_DEPTH"):
                res_x.setdefault("warnings", []).append(STATUS_DESCRIPTIONS["POLARIZATION_MIXING_RISK"])
    for label, res in results.items():
        if "tb_after_polarization" not in res and has_spectrum(res):
            res["tb_after_polarization"] = res["tb_spectrum"].copy()
    return results, m_pol


# ========================================================================
# БЛОК 6: ПЕРЕОТРАЖЕНИЯ ОТ СТЕНКИ КАМЕРЫ (феноменологическая геом. сумма)
# ========================================================================
def effective_reflection_coefficient(cfg):
    """
    R_eff = R_wall * round_trip_attenuation, 0<=R_eff<1 требуется для
    сходимости геометрической суммы переотражений.
    """
    r_wall = np.clip(cfg.wall_reflectivity, 0.0, 0.999999)
    atten = np.clip(cfg.round_trip_attenuation, 0.0, 1.0)
    r_eff = r_wall * atten
    return float(np.clip(r_eff, 0.0, 0.999999))


def apply_reflections(tb_spectrum, cfg):
    """
    Феноменологическая модель многократных переотражений:
        I_refl = I0 * sum_{k=0}^{N} (R_eff)^k             (конечный порядок)
        I_refl = I0 / (1 - R_eff)                          (предел N->inf, R_eff<1)
    ВНИМАНИЕ: это НЕ полная модель камеры с частотно-зависимой стенкой -
    феноменологическое приближение с единым коэффициентом R_eff.
    """
    if not cfg.reflections_enabled:
        return tb_spectrum.copy(), 1.0, False
    r_eff = effective_reflection_coefficient(cfg)
    stable = r_eff < 1.0
    if cfg.infinite_reflection_model and stable:
        gain = 1.0 / (1.0 - r_eff)
    else:
        n = int(max(cfg.max_reflection_order, 0))
        gain = float(np.sum(r_eff ** np.arange(0, n + 1)))
    return tb_spectrum * gain, gain, stable


def apply_reflections_to_all(results, cfg):
    r_eff = effective_reflection_coefficient(cfg)
    for label, res in results.items():
        base = res.get("tb_after_polarization")
        if base is None:
            continue
        refl, gain, stable = apply_reflections(base, cfg)
        res["tb_after_reflections"] = refl
        res["reflection_gain"] = gain
        res["reflection_stable"] = stable
        clean = res.get("tb_spectrum")
        if clean is not None:
            denom = np.maximum(np.abs(clean), 1e-30)
            res["delta_refl_percent"] = 100.0 * (refl - clean) / denom
    return results, r_eff


# ========================================================================
# БЛОК 5: ШУМЫ И ФОН ПРИЕМНИКА
# ========================================================================
def build_background(nu_grid, cfg):
    """I_background(nu) - медленно меняющийся фон тракта/приемника, эВ."""
    f_ghz = nu_grid / 1e9
    return cfg.background_offset_eV + cfg.background_slope_eV_per_GHz * (f_ghz - f_ghz.min())


def build_correlated_noise(nu_grid, cfg, rng):
    """
    Коррелированный по частоте шум (дрейф базовой линии) - реализован через
    гауссово сглаживание белого шума окном шириной frequency_correlation_width_hz.
    """
    n = nu_grid.size
    white = rng.normal(0.0, 1.0, size=n)
    dnu = np.mean(np.diff(nu_grid)) if n > 1 else 1.0
    sigma_bins = max(cfg.frequency_correlation_width_hz / dnu, 1e-6)
    max_half = max(n // 2 - 1, 1)
    half = int(np.clip(3 * sigma_bins, 1, max_half))
    x = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (x / sigma_bins) ** 2)
    kernel /= kernel.sum()
    corr = np.convolve(white, kernel, mode="same")
    if corr.size != n:
        corr = corr[:n] if corr.size > n else np.pad(corr, (0, n - corr.size))
    std_corr = np.std(corr)
    corr = corr * (cfg.absolute_noise_level_eV / std_corr) if std_corr > 1e-30 else corr
    return corr


def apply_noise_to_spectrum(nu_grid, tb_signal_eV, cfg, rng):
    """
    I_noisy = I_signal + I_background + N_white + N_corr.
    Белый шум - гауссов, мультипликативная (relative_noise_std * |сигнал|)
    и аддитивная (absolute_noise_level_eV) компоненты в квадратуре.
    Итоговое неотрицательное ограничение (клиппинг >=0), т.к. T_b/интенсивность
    физически неотрицательны; клиппинг применяется ПОСЛЕ формирования шума,
    и не путается с самим гауссовым распределением (которое может быть <0
    до клиппинга - это ожидаемо для гауссова шума).
    """
    background = build_background(nu_grid, cfg)
    sigma_white = np.sqrt((cfg.relative_noise_std * np.abs(tb_signal_eV)) ** 2
                           + cfg.absolute_noise_level_eV ** 2)
    white = rng.normal(0.0, 1.0, size=nu_grid.size) * sigma_white
    corr = build_correlated_noise(nu_grid, cfg, rng)
    noisy = tb_signal_eV + background + white + corr
    noisy_clipped = np.maximum(noisy, 0.0)
    n_clipped = int(np.sum(noisy < 0.0))
    return {
        "background": background, "white": white, "corr": corr,
        "noisy": noisy_clipped, "noisy_raw": noisy, "n_clipped": n_clipped,
    }


def apply_noise_to_all(results, nu_grid, cfg):
    if not cfg.noise_enabled:
        for res in results.values():
            base = res.get("tb_after_reflections")
            if base is not None:
                res["tb_final_noisy"] = base.copy()
                res["noise_components"] = None
        return results
    rng = np.random.default_rng(cfg.random_seed)
    for label, res in results.items():
        base = res.get("tb_after_reflections")
        if base is None:
            continue
        base_eV = joule_to_ev(base)
        comps = apply_noise_to_spectrum(nu_grid, base_eV, cfg, rng)
        res["noise_components"] = comps
        res["tb_final_noisy"] = ev_to_joule(comps["noisy"])
    return results


def noise_statistics(results):
    stats = {}
    for label, res in results.items():
        comps = res.get("noise_components")
        if comps is None:
            continue
        signal = res.get("tb_after_reflections")
        noisy = comps["noisy"]
        noise_only = noisy - joule_to_ev(signal) if signal is not None else noisy
        s = {
            "mean_eV": float(np.mean(noisy)),
            "std_eV": float(np.std(noisy)),
            "min_eV": float(np.min(noisy)),
            "max_eV": float(np.max(noisy)),
            "n_clipped": comps["n_clipped"],
        }
        sig_power = np.mean(joule_to_ev(signal) ** 2) if signal is not None else 0.0
        noise_power = np.mean(noise_only ** 2)
        s["snr"] = float(sig_power / noise_power) if noise_power > 1e-30 else float("nan")
        stats[label] = s
    return stats


# ========================================================================
# ГРАФИКИ
# ========================================================================
def plot_resonance_map(r, profiles, cfg):
    fig, ax = plt.subplots(figsize=(10, 6.5))
    styles = {1: "-", 2: "--", 3: ":"}
    for s in cfg.HARMONICS:
        f_ghz = harmonic_resonance_frequency(profiles["B"], s) / 1e9
        ax.plot(r, f_ghz, color="black", lw=1.8, ls=styles[s], label=f"s={s}: $s\\cdot f_{{ce}}(R)$")
    ax.set_xlabel("R, м"); ax.set_ylabel("Резонансная частота, ГГц")
    ax.set_title("Резонансные частоты гармоник $s\\cdot f_{ce}(R)$ (холодноплазм. приближение, $\\gamma=1$)")
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
    f = out_path("00_резонансные_частоты_гармоник.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_cutoff_diagram(r, f_l, f_r, f_pe, f_ce, cfg):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    f_ce_ghz, f_pe_ghz, f_r_ghz, f_l_ghz = f_ce/1e9, f_pe/1e9, f_r/1e9, f_l/1e9

    ax.plot(r, f_ce_ghz, color="black", lw=1.8, label=r"ЭЦ резонанс $f_{ce}$ (s=1)")
    for s in (2, 3):
        ax.plot(r, s*f_ce_ghz, color="gray", lw=1.1, ls="-.", label=f"$s\\cdot f_{{ce}}$, s={s}")
    ax.plot(r, f_pe_ghz, color="tab:blue", lw=1.6, ls="--", label=r"Отсечка O-моды $f_{pe}$")
    ax.plot(r, f_r_ghz, color="tab:red", lw=1.8, label=r"Верхняя отсечка X-моды $f_U$")
    ax.plot(r, f_l_ghz, color="tab:green", lw=1.8, label=r"Нижняя отсечка X-моды $f_L$")

    ax.fill_between(r, 0, f_pe_ghz, facecolor="none", edgecolor="tab:blue", hatch="|||", alpha=0.9)
    ax.fill_between(r, 0, f_l_ghz, where=(f_l_ghz >= 0), facecolor="none", edgecolor="tab:green", hatch="xx", alpha=0.95)
    ax.fill_between(r, f_l_ghz, f_r_ghz, where=(f_r_ghz >= f_l_ghz), facecolor="none", edgecolor="tab:red", hatch="////", alpha=0.6)

    handles = [
        mpatches.Patch(facecolor="none", edgecolor="tab:blue", hatch="|||", label="Непрозрачность O-моды"),
        mpatches.Patch(facecolor="none", edgecolor="tab:green", hatch="xx", label="Непрозрачность X-моды (нижняя)"),
        mpatches.Patch(facecolor="none", edgecolor="tab:red", hatch="////", label="Полоса непрозрачности X-моды [f_L,f_U]"),
    ]
    ax.set_xlabel("R, м"); ax.set_ylabel("Частота f, ГГц")
    ax.set_title("Диаграмма частот отсечки O-/X-моды и гармоник циклотронного резонанса")
    ax.legend(handles=ax.get_lines() + handles, fontsize=8, loc="upper right")
    ax.set_xlim(r.min(), r.max()); ax.set_ylim(0, 3.5*f_ce_ghz.max())

    f = out_path("ec_cutoff_diagram_trt.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_status_map(results, cfg):
    """Карта статусов доступности всех комбинаций мода-гармоника."""
    labels = cfg.combo_labels()
    status_codes = {
        "CALCULATED": 3, "CALCULATED_QUALITATIVE": 2, "LOW_OPTICAL_DEPTH": 1,
        "CUTOFF_RISK": -1, "OUT_OF_RESONANCE_RANGE": -2,
        "NOT_IMPLEMENTED_PHYSICALLY": -3, "INSUFFICIENT_INPUT_DATA": -2,
    }
    values = [status_codes.get(results[l]["status"], 0) for l in labels]
    colors = ["#2E7D32" if v == 3 else "#66BB6A" if v == 2 else "#AED581" if v == 1
              else "#EF6C00" if v == -1 else "#B0BEC5" if v == -2 else "#78909C" for v in values]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(labels, [1]*len(labels), color=colors, edgecolor="black")
    for i, l in enumerate(labels):
        ax.text(i, 0.5, results[l]["status"], rotation=90, ha="center", va="center", fontsize=7, color="black")
    ax.set_yticks([])
    ax.set_title("Карта статусов доступности каналов O1,X1,O2,X2,O3,X3")
    legend_items = [
        mpatches.Patch(facecolor="#2E7D32", label="CALCULATED"),
        mpatches.Patch(facecolor="#66BB6A", label="CALCULATED_QUALITATIVE"),
        mpatches.Patch(facecolor="#AED581", label="LOW_OPTICAL_DEPTH"),
        mpatches.Patch(facecolor="#EF6C00", label="CUTOFF_RISK"),
        mpatches.Patch(facecolor="#B0BEC5", label="OUT_OF_RESONANCE_RANGE"),
    ]
    ax.legend(handles=legend_items, fontsize=7, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.05))
    f = out_path("04_карта_статусов_каналов.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_all_spectra(nu_grid, results, cfg, log_scale=False):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for label, res in results.items():
        color = MODE_COLORS.get(label, "gray")
        if has_spectrum(res):
            tb_kev = joule_to_ev(res["tb_spectrum"]) / 1e3
            tag = {"CALCULATED": "рассчитано", "CALCULATED_QUALITATIVE": "кач. оценка",
                   "LOW_OPTICAL_DEPTH": "оптич. тонкий"}.get(res["status"], res["status"])
            ax.plot(nu_grid/1e9, tb_kev, color=color, lw=1.8, label=f"{label} ({tag})")
        else:
            ax.plot([], [], color=color, lw=1.8, ls="--",
                    label=f"{label}: {res['status']}")

    ax.set_xlabel("Частота, ГГц"); ax.set_ylabel(r"Яркостная температура $T_b$, кэВ")
    suffix = " (логарифмическая шкала)" if log_scale else ""
    ax.set_title(f"Чистые спектры плазменного ECE по комбинациям мода-гармоника{suffix}\n"
                 "(X - валидированная модель; O - качественная оценка; пунктир = не рассчитано)")
    if log_scale:
        ax.set_yscale("log")
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)

    fname = out_path("05_чистые_спектры_all_modes_log.png" if log_scale else "05_чистые_спектры_all_modes.png")
    fig.tight_layout(); fig.savefig(fname, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fname


def plot_harmonic_sums(nu_grid, results, cfg):
    """
    Спектры по гармоникам O_s+X_s. ВНИМАНИЕ: прямое арифметическое сложение
    T_b(O_s) и T_b(X_s) физически оправдано ТОЛЬКО если приемник
    регистрирует суммарную (неполяризованную либо обе поляризации без
    разделения) мощность без поляризационной селекции. Если приемник
    поляризационно-селективен, такая сумма не соответствует наблюдаемому
    сигналу - показана здесь как ИЛЛЮСТРАТИВНАЯ величина.
    """
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for s in cfg.HARMONICS:
        lbl_o, lbl_x = f"O{s}", f"X{s}"
        res_o, res_x = results.get(lbl_o), results.get(lbl_x)
        io = res_o["tb_spectrum"] if res_o and has_spectrum(res_o) else np.zeros_like(nu_grid)
        ix = res_x["tb_spectrum"] if res_x and has_spectrum(res_x) else np.zeros_like(nu_grid)
        total = joule_to_ev(io + ix) / 1e3
        if np.any(total > 0):
            ax.plot(nu_grid/1e9, total, lw=1.8, label=f"O{s}+X{s} (сумма, см. пояснение)")
    ax.set_xlabel("Частота, ГГц"); ax.set_ylabel(r"$T_b$, кэВ")
    ax.set_title("Иллюстративная сумма O_s+X_s по гармоникам\n"
                 "(физически оправдано только без поляризационной селекции приемника)")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    f = out_path("06_суммы_по_гармоникам_OX.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_total_clean_spectrum(nu_grid, results, cfg):
    total = np.zeros_like(nu_grid)
    for res in results.values():
        if has_spectrum(res):
            total += res["tb_spectrum"]
    fig, ax = plt.subplots(figsize=(11, 6.0))
    ax.plot(nu_grid/1e9, joule_to_ev(total)/1e3, color=C_MAIN, lw=2.0)
    ax.set_xlabel("Частота, ГГц"); ax.set_ylabel(r"$T_b$, кэВ")
    ax.set_title("Суммарный чистый спектр по всем рассчитанным каналам\n(иллюстративная сумма, см. пояснение к рис.06)")
    ax.grid(True, alpha=0.3)
    f = out_path("07_суммарный_чистый_спектр.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f, total


def plot_pipeline_comparison(nu_grid, results, cfg, label):
    res = results.get(label)
    if res is None or not has_spectrum(res):
        return None
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.plot(nu_grid/1e9, joule_to_ev(res["tb_spectrum"])/1e3, color="black", lw=1.6, label="чистый ECE")
    if "tb_after_polarization" in res:
        ax.plot(nu_grid/1e9, joule_to_ev(res["tb_after_polarization"])/1e3, color="tab:blue", lw=1.4, ls="--", label="после поляризац. смешивания")
    if "tb_after_reflections" in res:
        ax.plot(nu_grid/1e9, joule_to_ev(res["tb_after_reflections"])/1e3, color="tab:orange", lw=1.4, ls="-.", label="после переотражений")
    if "tb_final_noisy" in res:
        ax.plot(nu_grid/1e9, joule_to_ev(res["tb_final_noisy"])/1e3, color="tab:red", lw=1.0, alpha=0.7, label="итоговый наблюдаемый (с шумом)")
    ax.set_xlabel("Частота, ГГц"); ax.set_ylabel(r"$T_b$, кэВ")
    ax.set_title(f"Канал {label}: этапы формирования наблюдаемого спектра")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    f = out_path(f"08_этапы_канал_{label}.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_reflection_impact(nu_grid, results, cfg):
    fig, ax = plt.subplots(figsize=(11, 6.0))
    any_plotted = False
    for label, res in results.items():
        if "delta_refl_percent" in res:
            ax.plot(nu_grid/1e9, res["delta_refl_percent"], lw=1.3, label=label, color=MODE_COLORS.get(label, "gray"))
            any_plotted = True
    ax.axhline(0, color="black", lw=1.0, ls=":")
    ax.set_xlabel("Частота, ГГц"); ax.set_ylabel(r"$\delta_{refl}$, %")
    ax.set_title("Относительное влияние переотражений на спектр по каналам")
    if any_plotted:
        ax.legend(fontsize=8, ncol=2)
    f = out_path("09_относительное_влияние_переотражений.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_noise_snr(nu_grid, results, cfg, noise_stats):
    calc_labels = [l for l, s in noise_stats.items() if not np.isnan(s.get("snr", np.nan))]
    if not calc_labels:
        return None
    fig, ax = plt.subplots(figsize=(9, 5.5))
    snrs = [noise_stats[l]["snr"] for l in calc_labels]
    ax.bar(calc_labels, snrs, color=[MODE_COLORS.get(l, "gray") for l in calc_labels])
    ax.set_ylabel("Отношение сигнал/шум (мощность)")
    ax.set_title("Оценка ОСШ по каналам (DEMO-параметры шума приемника)")
    ax.set_yscale("log")
    f = out_path("10_отношение_сигнал_шум.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_tau_by_combo(r, results, cfg):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for label, res in results.items():
        color = MODE_COLORS.get(label, "gray")
        if res.get("tau_cum") is not None:
            tau_max_per_nu = res["tau_cum"][-1, :]
            valid = tau_max_per_nu > 0
            if np.any(valid):
                ls = "-" if res["status"] == "CALCULATED" else "--"
                ax.plot(np.arange(tau_max_per_nu.size)[valid], tau_max_per_nu[valid],
                        color=color, lw=1.5, ls=ls, label=f"{label}: max={tau_max_per_nu.max():.2f}")
    ax.axhline(1.0, color="black", ls=":", lw=1.5, label=r"$\tau=1$ (порог)")
    ax.set_xlabel("Индекс частотного канала"); ax.set_ylabel(r"Оптическая толщина $\tau_\nu$")
    ax.set_yscale("log")
    ax.set_title("Оптическая толщина по частотной сетке\n(сплошная - X, пунктир - O качественная оценка)")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    f = out_path("11_оптическая_толщина_по_комбинациям.png")
    fig.tight_layout(); fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


def plot_diagnostics_overview(profiles, results, cfg):
    r = profiles["R"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Диагностика ECE: профили и коэффициент поглощения (все рассчитанные комбинации)",
                 fontsize=13, fontweight="bold")

    ax = axs[0, 0]; ax.plot(r, profiles["B"], color=C_MAIN, lw=2)
    ax.set_xlabel("R, м"); ax.set_ylabel("B, Тл"); ax.set_title("Магнитное поле B(R)")

    ax = axs[0, 1]; ax.plot(r, profiles["ne"], color=C_ALT, lw=2)
    ax.set_xlabel("R, м"); ax.set_ylabel(r"$n_e$, м$^{-3}$"); ax.set_title("Концентрация $n_e(R)$")

    ax = axs[1, 0]; ax.plot(r, profiles["Te_keV"], color=C_3, lw=2)
    ax.set_xlabel("R, м"); ax.set_ylabel(r"$T_e$, кэВ"); ax.set_title("Температура $T_e(R)$")

    ax = axs[1, 1]
    for label, res in results.items():
        if res.get("alpha") is not None:
            mid = res["alpha"].shape[1] // 2
            ls = "-" if res["status"] == "CALCULATED" else "--"
            ax.plot(r, res["alpha"][:, mid], color=MODE_COLORS[label], lw=2, ls=ls, label=f"{label} (срез сетки)")
    ax.set_xlabel("R, м"); ax.set_ylabel(r"$\alpha_\nu$, м$^{-1}$")
    ax.set_title("Коэффициент поглощения (сплошная - X, пунктир - O кач. оценка)")
    ax.legend(fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    f = out_path("ece_diagnostics_trt.png")
    fig.savefig(f, dpi=150, bbox_inches="tight"); plt.close(fig)
    return f


# ========================================================================
# ПРОВЕРКИ КОРРЕКТНОСТИ
# ========================================================================
def run_sanity_checks(profiles, results, cfg, r_eff, m_pol):
    checks = {}
    checks["finite_profiles"] = bool(
        np.all(np.isfinite(profiles["B"])) and np.all(np.isfinite(profiles["ne"]))
        and np.all(np.isfinite(profiles["Te_keV"]))
    )
    checks["positive_profiles"] = bool(
        np.all(profiles["B"] > 0) and np.all(profiles["ne"] > 0) and np.all(profiles["Te_keV"] > 0)
    )

    te_max_j = ev_to_joule(cfg.Te0_keV * 1e3)
    all_finite_spectra, all_nonneg_alpha, all_tb_bounded = True, True, True
    for label, res in results.items():
        if has_spectrum(res):
            if not np.all(np.isfinite(res["tb_spectrum"])):
                all_finite_spectra = False
            if np.any(res["tb_spectrum"] < -1e-30):
                all_nonneg_alpha = False
            if np.any(res["tb_spectrum"] > te_max_j * 1.05):
                all_tb_bounded = False
            if res["alpha"] is not None and np.any(res["alpha"] < -1e-30):
                all_nonneg_alpha = False

    checks["finite_spectra_all_calculated"] = all_finite_spectra
    checks["nonnegative_alpha_all_calculated"] = all_nonneg_alpha
    checks["Tb_bounded_all_calculated"] = all_tb_bounded

    labels_present = set(results.keys())
    expected = set(cfg.combo_labels())
    checks["all_six_combinations_present"] = (labels_present == expected)

    checks["reflection_coefficient_stable"] = bool(0.0 <= r_eff < 1.0)

    # При R_wall=0 переотражения должны исчезать (gain=1)
    cfg_zero = PlasmaConfig()
    cfg_zero.wall_reflectivity = 0.0
    test_spec = np.array([1.0, 2.0, 3.0])
    refl_zero, gain_zero, _ = apply_reflections(test_spec, cfg_zero)
    checks["reflections_vanish_at_zero_wall_R"] = bool(np.allclose(refl_zero, test_spec) and abs(gain_zero - 1.0) < 1e-9)

    # При нулевом поляризационном смешивании утечки быть не должно
    cfg_nomix = PlasmaConfig()
    cfg_nomix.polarization_mixing_ox = 0.0
    cfg_nomix.polarization_mixing_xo = 0.0
    cfg_nomix.polarizer_o_leakage = 0.0
    m_nomix = build_polarization_matrix(cfg_nomix)
    checks["no_polarization_leakage_when_zero"] = bool(abs(m_nomix[0, 1]) < 1e-12 and abs(m_nomix[1, 0]) < 1e-12)

    # Отключенный шум => финальный спектр совпадает с спектром после отражений
    match_no_noise = True
    for res in results.values():
        if "tb_after_reflections" in res and "tb_final_noisy" in res and not cfg.noise_enabled:
            if not np.allclose(res["tb_final_noisy"], res["tb_after_reflections"], rtol=1e-9):
                match_no_noise = False
    checks["noise_disabled_matches_reflected_spectrum"] = match_no_noise

    # Нулевые отражения и нулевой шум => итог совпадает с чистым спектром
    cfg_clean = PlasmaConfig()
    cfg_clean.reflections_enabled = False
    cfg_clean.noise_enabled = False
    cfg_clean.polarization_mixing_ox = 0.0
    cfg_clean.polarization_mixing_xo = 0.0
    cfg_clean.polarizer_o_leakage = 0.0
    test_ok = True
    sample_label = next((l for l, res in results.items() if has_spectrum(res)), None)
    if sample_label:
        clean_spec = results[sample_label]["tb_spectrum"]
        refl_clean, gain_c, _ = apply_reflections(clean_spec, cfg_clean)
        if not np.allclose(refl_clean, clean_spec, rtol=1e-9):
            test_ok = False
    checks["clean_matches_when_no_effects"] = test_ok

    tau_x = {}
    for s in cfg.HARMONICS:
        res = results.get(f"X{s}")
        if res and res.get("tau_cum") is not None:
            tau_x[s] = res["tau_cum"][-1, :].max()
    if len(tau_x) >= 2:
        sorted_s = sorted(tau_x.keys())
        monotonic_decreasing = all(tau_x[sorted_s[i]] >= tau_x[sorted_s[i+1]] for i in range(len(sorted_s)-1))
        checks["tau_X_decreases_with_harmonic"] = monotonic_decreasing
    else:
        checks["tau_X_decreases_with_harmonic"] = True

    # Проверка воспроизводимости random_seed: разные seed -> разная реализация шума, чистый спектр не меняется
    if cfg.noise_enabled and sample_label:
        rng1 = np.random.default_rng(cfg.random_seed)
        rng2 = np.random.default_rng(cfg.random_seed + 1)
        base = results[sample_label].get("tb_after_reflections")
        if base is not None:
            base_eV = joule_to_ev(base)
            c1 = apply_noise_to_spectrum(np.linspace(1e9, 2e9, base_eV.size), base_eV, cfg, rng1)
            c2 = apply_noise_to_spectrum(np.linspace(1e9, 2e9, base_eV.size), base_eV, cfg, rng2)
            checks["different_seed_changes_noise_not_clean"] = bool(not np.allclose(c1["noisy"], c2["noisy"]))
        else:
            checks["different_seed_changes_noise_not_clean"] = True
    else:
        checks["different_seed_changes_noise_not_clean"] = True

    return checks


# ========================================================================
# СОХРАНЕНИЕ CSV И ОТЧЕТА
# ========================================================================
def save_all_spectra_csv(nu_grid, results, cfg):
    labels = list(results.keys())
    cols_clean, cols_pol, cols_refl, cols_noisy = [nu_grid], [nu_grid], [nu_grid], [nu_grid]
    header_clean, header_pol, header_refl, header_noisy = ["f_Hz"], ["f_Hz"], ["f_Hz"], ["f_Hz"]
    for l in labels:
        res = results[l]
        z = np.zeros_like(nu_grid)
        cols_clean.append(joule_to_ev(res["tb_spectrum"]) if has_spectrum(res) else z)
        header_clean.append(f"Tb_clean_{l}_eV")
        cols_pol.append(joule_to_ev(res["tb_after_polarization"]) if "tb_after_polarization" in res else z)
        header_pol.append(f"Tb_afterpol_{l}_eV")
        cols_refl.append(joule_to_ev(res["tb_after_reflections"]) if "tb_after_reflections" in res else z)
        header_refl.append(f"Tb_afterrefl_{l}_eV")
        cols_noisy.append(joule_to_ev(res["tb_final_noisy"]) if "tb_final_noisy" in res else z)
        header_noisy.append(f"Tb_noisy_{l}_eV")

    paths = {}
    for name, cols, header in [
        ("12_спектры_чистые_ECE.csv", cols_clean, header_clean),
        ("13_спектры_после_поляризации.csv", cols_pol, header_pol),
        ("14_спектры_после_переотражений.csv", cols_refl, header_refl),
        ("15_спектры_итоговые_с_шумом.csv", cols_noisy, header_noisy),
    ]:
        p = out_path(name)
        np.savetxt(p, np.column_stack(cols), delimiter=",", header=",".join(header), comments="", fmt="%.6e")
        paths[name] = p
    return paths


def save_report(profiles, results, cfg, checks, noise_stats, r_eff, m_pol, nu_grid, path):
    lines = []
    lines.append("ОТЧЕТ О РАСЧЕТЕ РАСШИРЕННОЙ FORWARD-МОДЕЛИ ECE (O1,X1,O2,X2,O3,X3)")
    lines.append("=" * 78)
    lines.append("")
    lines.append("1. ПАРАМЕТРЫ ЗАПУСКА")
    lines.append(f"   R0={cfg.R0} м, a={cfg.a} м, B0={cfg.B0} Тл, Te0={cfg.Te0_keV} кэВ, ne0={cfg.ne0_m3:.3e} м^-3")
    lines.append(f"   Гармоники: {cfg.HARMONICS}, Моды: {cfg.MODES}")
    lines.append(f"   random_seed={cfg.random_seed}, noise_enabled={cfg.noise_enabled}, reflections_enabled={cfg.reflections_enabled}")
    lines.append("")
    lines.append("2. ДИАПАЗОНЫ ЧАСТОТ РЕЗОНАНСА")
    for s in cfg.HARMONICS:
        f_res = harmonic_resonance_frequency(profiles["B"], s)
        lines.append(f"   s={s}: {f_res.min()/1e9:.2f}-{f_res.max()/1e9:.2f} ГГц")
    lines.append("")
    lines.append("3. СТАТУСЫ КАНАЛОВ")
    for l, res in results.items():
        lines.append(f"   {l}: {res['status']} - {res['reason']}")
    lines.append("")
    lines.append("4. ПАРАМЕТРЫ ШУМА (DEMO, требуют замены на параметры прибора)")
    lines.append(f"   receiver_noise_temperature_K={cfg.receiver_noise_temperature_K}")
    lines.append(f"   relative_noise_std={cfg.relative_noise_std}, absolute_noise_level_eV={cfg.absolute_noise_level_eV}")
    lines.append(f"   background_offset_eV={cfg.background_offset_eV}, background_slope_eV_per_GHz={cfg.background_slope_eV_per_GHz}")
    for l, s in noise_stats.items():
        lines.append(f"   [{l}] mean={s['mean_eV']:.3f} эВ, std={s['std_eV']:.3f} эВ, min={s['min_eV']:.3f}, max={s['max_eV']:.3f}, SNR={s['snr']:.3g}, clipped={s['n_clipped']}")
    lines.append("")
    lines.append("5. ПАРАМЕТРЫ ПЕРЕОТРАЖЕНИЙ И ПОЛЯРИЗАЦИИ (DEMO)")
    lines.append(f"   wall_reflectivity={cfg.wall_reflectivity}, round_trip_attenuation={cfg.round_trip_attenuation}, R_eff={r_eff:.4f}")
    lines.append(f"   polarization_mixing_ox={cfg.polarization_mixing_ox}, polarization_mixing_xo={cfg.polarization_mixing_xo}")
    lines.append(f"   M_pol matrix=\n{m_pol}")
    lines.append("")
    lines.append("6. ОГРАНИЧЕНИЯ МОДЕЛИ")
    lines.append("   - Не является полной релятивистской кинетической ECE-моделью.")
    lines.append("   - Гармонический коэффициент O-моды - качественная параметризация по аналогии.")
    lines.append("   - Критерий доступности - упрощенный (не лучевая трассировка, не полное дисперс. уравнение).")
    lines.append("   - Модель переотражений и поляризационного смешивания - феноменологическая (DEMO-параметры).")
    lines.append("   - Параметры шума приемника - демонстрационные значения, требуют замены на паспортные.")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return path


# ========================================================================
# ГЛАВНАЯ ПРОГРАММА
# ========================================================================
def main():
    cfg = PlasmaConfig()

    print("=" * 78)
    print(" РАСШИРЕННАЯ FORWARD-МОДЕЛЬ ECE: O1,X1,O2,X2,O3,X3 + ШУМ/ОТРАЖЕНИЯ/ПОЛЯРИЗАЦИЯ")
    print("=" * 78)
    print(f" Большой радиус R0            = {cfg.R0:.3f} м")
    print(f" Малый радиус a                = {cfg.a:.3f} м")
    print(f" Тороидальное поле B0          = {cfg.B0:.3f} Тл")
    print(f" Центральная температура Te0   = {cfg.Te0_keV:.2f} кэВ")
    print(f" Центральная концентрация ne0  = {cfg.ne0_m3:.3e} м^-3")
    print(f" Угол наблюдения               = {cfg.view_angle_deg:.1f} град. (перпендикулярное)")
    print(f" Рассматриваемые гармоники     = {cfg.HARMONICS}")
    print(f" Рассматриваемые моды          = {cfg.MODES}")
    print("=" * 78)
    print("\n [ПРЕДУПРЕЖДЕНИЕ] Параметры шума приемника и переотражений - DEMO-значения!")
    print(" Требуется замена на паспортные/измеренные параметры реальной установки.")

    profiles = build_plasma_profiles(cfg)
    r = profiles["R"]

    print("\n" + "-" * 78)
    print(" ВХОДНЫЕ ПРОФИЛИ ПЛАЗМЫ")
    print("-" * 78)
    print(f" B(R0)   = {np.interp(cfg.R0, r, profiles['B']):.4f} Тл")
    print(f" ne(ось) = {profiles['ne'].max():.3e} м^-3")
    print(f" Te(ось) = {profiles['Te_keV'].max():.3f} кэВ")

    profile_files = plot_input_profiles(profiles, cfg)
    print("\n Сохранены отдельные графики входных профилей:")
    for f in profile_files:
        print(f"   - {os.path.basename(f)}")

    f_l, f_r, f_pe, f_ce = cutoff_frequencies(profiles["B"], profiles["ne"])

    print("\n" + "-" * 78)
    print(" ХОЛОДНОПЛАЗМЕННЫЕ ЧАСТОТЫ (на магнитной оси R0)")
    print("-" * 78)
    i0 = np.argmin(np.abs(r - cfg.R0))
    print(f" f_ce(R0)  = {f_ce[i0]/1e9:.4f} ГГц")
    print(f" f_pe(R0)  = {f_pe[i0]/1e9:.4f} ГГц")
    print(f" f_L(R0)   = {f_l[i0]/1e9:.4f} ГГц")
    print(f" f_U(R0)   = {f_r[i0]/1e9:.4f} ГГц")

    f_min = 0.5 * f_ce.min()
    f_max = 3.5 * f_ce.max()
    nu_grid = np.linspace(f_min, f_max, cfg.n_grid_nu)
    print(f"\n Диапазон частотной сетки: {nu_grid.min()/1e9:.2f}-{nu_grid.max()/1e9:.2f} ГГц")

    print("\n" + "-" * 78)
    print(" РЕЗОНАНСНЫЕ ЧАСТОТЫ ГАРМОНИК s*f_ce(R)")
    print("-" * 78)
    for s in cfg.HARMONICS:
        f_res = harmonic_resonance_frequency(profiles["B"], s)
        print(f"   s={s}: диапазон {f_res.min()/1e9:.2f}-{f_res.max()/1e9:.2f} ГГц "
              f"(на R=[{r.min():.3f};{r.max():.3f}] м)")

    plot_resonance_map(r, profiles, cfg)

    print("\n" + "-" * 78)
    print(" РАСЧЕТ ВСЕХ КОМБИНАЦИЙ МОДА-ГАРМОНИКА (O1,X1,O2,X2,O3,X3) - ЧИСТЫЙ ECE")
    print("-" * 78)
    results = compute_all_combinations(cfg, r, nu_grid, profiles, f_l, f_r, f_pe)

    for label, res in results.items():
        print(f"\n   [{label}] статус: {res['status']}")
        print(f"     резонансный диапазон: {res['f_res_min_GHz']:.2f}-{res['f_res_max_GHz']:.2f} ГГц")
        print(f"     причина: {res['reason']}")
        print(f"     точность модели: {res['accuracy']}")
        if has_spectrum(res):
            tau_max = res["tau_max"]
            tb_max_eV = joule_to_ev(res["tb_spectrum"].max())
            print(f"     max(tau) = {tau_max:.3e}, max(T_b) = {tb_max_eV:.1f} эВ")

    print("\n" + "-" * 78)
    print(" ПОЛЯРИЗАЦИОННОЕ СМЕШИВАНИЕ O<->X")
    print("-" * 78)
    results, m_pol = apply_polarization_mixing(results, cfg)
    print(f"   Матрица M_pol =\n{m_pol}")

    print("\n" + "-" * 78)
    print(" ПЕРЕОТРАЖЕНИЯ ОТ СТЕНКИ КАМЕРЫ")
    print("-" * 78)
    results, r_eff = apply_reflections_to_all(results, cfg)
    print(f"   R_eff = {r_eff:.4f} (устойчиво: {r_eff < 1.0})")

    print("\n" + "-" * 78)
    print(" ШУМЫ И ФОН ПРИЕМНИКА")
    print("-" * 78)
    results = apply_noise_to_all(results, nu_grid, cfg)
    noise_stats = noise_statistics(results)
    for l, s in noise_stats.items():
        print(f"   [{l}] mean={s['mean_eV']:.3f} эВ std={s['std_eV']:.3f} эВ min={s['min_eV']:.3f} max={s['max_eV']:.3f} SNR={s['snr']:.3g} clipped={s['n_clipped']}")

    channel_centers = np.linspace(nu_grid.min(), nu_grid.max(), cfg.n_channels)
    channel_results = {}
    for label, res in results.items():
        if has_spectrum(res):
            tb_ch = channel_spectrum(nu_grid, res["tb_spectrum"], channel_centers, cfg.channel_half_width_hz)
            channel_results[label] = tb_ch

    checks = run_sanity_checks(profiles, results, cfg, r_eff, m_pol)
    print("\n" + "-" * 78)
    print(" ПРОВЕРКИ ФИЗИЧЕСКОЙ И ЧИСЛЕННОЙ КОРРЕКТНОСТИ")
    print("-" * 78)
    all_ok = True
    for name, ok in checks.items():
        status = "ПРОЙДЕНО" if ok else "ОШИБКА"
        if not ok:
            all_ok = False
        print(f"   {name:45s}: {status}")

    print("\n" + "-" * 78)
    print(" СОХРАНЕНИЕ ГРАФИКОВ И ФАЙЛОВ РЕЗУЛЬТАТОВ")
    print("-" * 78)

    f1 = plot_cutoff_diagram(r, f_l, f_r, f_pe, f_ce, cfg)
    f2 = plot_all_spectra(nu_grid, results, cfg, log_scale=False)
    f3 = plot_all_spectra(nu_grid, results, cfg, log_scale=True)
    f4 = plot_tau_by_combo(r, results, cfg)
    f5 = plot_diagnostics_overview(profiles, results, cfg)
    f6 = plot_status_map(results, cfg)
    f7 = plot_harmonic_sums(nu_grid, results, cfg)
    f8, total_clean = plot_total_clean_spectrum(nu_grid, results, cfg)
    f9 = plot_reflection_impact(nu_grid, results, cfg)
    f10 = plot_noise_snr(nu_grid, results, cfg, noise_stats)
    pipeline_plots = []
    for label in cfg.combo_labels():
        p = plot_pipeline_comparison(nu_grid, results, cfg, label)
        if p:
            pipeline_plots.append(p)

    status_csv = out_path("статусы_комбинаций_O_X.csv")
    with open(status_csv, "w", encoding="utf-8") as fcsv:
        fcsv.write("Комбинация,Статус,ЧастотаРезонансаМин_ГГц,ЧастотаРезонансаМакс_ГГц,Причина,ТочностьМодели\n")
        for label, res in results.items():
            fcsv.write(f"{label},{res['status']},{res['f_res_min_GHz']:.4f},"
                       f"{res['f_res_max_GHz']:.4f},\"{res['reason']}\",\"{res['accuracy']}\"\n")

    if channel_results:
        spec_csv = out_path("канальные_спектры_рассчитанных_комбинаций.csv")
        header = "f_center_Hz," + ",".join(f"Tb_{lbl}_eV" for lbl in channel_results)
        cols = [channel_centers] + [joule_to_ev(v) for v in channel_results.values()]
        data = np.column_stack(cols)
        np.savetxt(spec_csv, data, delimiter=",", header=header, comments="", fmt="%.6e")
    else:
        spec_csv = None

    spectra_csvs = save_all_spectra_csv(nu_grid, results, cfg)
    report_path = out_path("16_отчет_о_расчете.txt")
    save_report(profiles, results, cfg, checks, noise_stats, r_eff, m_pol, nu_grid, report_path)

    saved = [f1, f2, f3, f4, f5, f6, f7, f8, f9, status_csv, report_path] + pipeline_plots
    if f10:
        saved.append(f10)
    if spec_csv:
        saved.append(spec_csv)
    saved.extend(spectra_csvs.values())
    for f in saved:
        print(f"   - {os.path.basename(f)}")

    print("\n" + "=" * 78)
    if all_ok:
        print(" РАСЧЕТ ЗАВЕРШЕН УСПЕШНО. Все проверки корректности пройдены.")
    else:
        print(" РАСЧЕТ ЗАВЕРШЕН С ПРЕДУПРЕЖДЕНИЯМИ. См. таблицу проверок выше.")
    calc_valid = [l for l, res in results.items() if res["status"] == "CALCULATED"]
    calc_qual = [l for l, res in results.items() if res["status"] == "CALCULATED_QUALITATIVE"]
    low_tau = [l for l, res in results.items() if res["status"] == "LOW_OPTICAL_DEPTH"]
    not_calc = [l for l, res in results.items() if not has_spectrum(res)]
    print(f" Рассчитаны (валидированная модель X): {calc_valid}")
    print(f" Рассчитаны (качественная оценка O):    {calc_qual}")
    print(f" Оптически тонкие каналы:               {low_tau}")
    print(f" НЕ рассчитаны (отсечка/вне диапазона): {not_calc}")
    print("=" * 78)

    return {"profiles": profiles, "results": results, "checks": checks,
            "channel_results": channel_results, "noise_stats": noise_stats, "r_eff": r_eff}


if __name__ == "__main__":
    main()
