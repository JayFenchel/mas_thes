#include <math.h>
#include <stdio.h>

float * func(float *x, float *u, float delta_t ) {

	// Physikalische Eigenschaften
	// Gemessen
	const float m = 2.;
	// [kg] Gewicht                      *nicht ueberprueft*
	const float r = 0.0335;
	// [m] Radradius                     *nicht ueberprueft*
	const float g = 9.80665;
	// [m/s^2] Erdbeschleunigung
	const float l = 0.2570;
	// [m] Gesamtlaenge des Radstandes   *nicht ueberprueft*
	const float rho_l = 1.2041;
	// [kg/m^3] Luftdichte (Meeresspiegel, 20 Grad Celsius)
	// Aus Ueberlegung geschaetzt
	const float m_distribution = 0.5;
	// Verhaeltnis lv/l (lv - Abstand Vorderachse zum Fahrzeugschwerpunkt)
	const float T = 0.2;
	// [s] zeitliche Verzoegerung des Inputs
	float J = (1. / 12. * m * (0.38 * 0.38 + 0.185 * 0.185));
	// [kg*m^2] Traegheitsmoment eines Qauders mit den Kantenlaengen a=0.38m, b=0.185m
	float c_wA = 1.15 * 0.025944; // [-]*[m^2] Stroemungswiderstandskoeffizient mal Flaeche

	// Aus Messdaten gschaetzt
	const float m_gamma = 0.5061;
	const float xi = 84.5;
	const float f_r = 10.6577;
	const float c_value = 24.4686;
	const float c_alpha_distribution = 0.3844;

	float xi_v = xi;
	float xi_h = xi;
	float c_h_alpha = c_value;

	float c_v_alpha = c_h_alpha * c_alpha_distribution;

	// Zusammengesetzte Parameter für Zustandsgleichungen
	float lv = l * m_distribution; // [m] Laenge vordere Haelfte
	float lh = l - lv; // [m] Laenge hintere Haelfte

	float p1 = c_h_alpha * lh / J;
	float p2 = c_v_alpha * lv / J;
	float p3 = m_gamma * lv / J / r / xi_v;
	float p4 = c_v_alpha / m;
	float p5 = m_gamma / m / r / xi_v;
	float p6 = c_h_alpha / m;
	float p7 = (1. - m_gamma) / r / xi_h / m;
	float p8 = f_r / m;
	float p9 = ((c_h_alpha * (m / 2.) * (m / 2.) + c_v_alpha * (m / 2.) * (m / 2.)) / m / c_v_alpha / c_h_alpha);
	float p10 = c_wA * rho_l / 2. / m;
	float p11 = T;

//// Meine Modellgleichungen aus Matlab
	float dx0 = (x[5] * cos(x[2] + x[4])) * delta_t + x[0];
	float dx1 = (x[5] * sin(x[2] + x[4])) * delta_t + x[1];
	float dx2 = (x[3]) * delta_t + x[2];

	float dx3 = (-p1 * atan(lh * x[3] / x[5]) + p2 * x[6] * cos(x[6]) - p2 * atan(lv * x[3] / x[5]) * cos(x[6]) + p3 * u[1] * x[6]) * delta_t + x[3];
	float dx4 = (-x[3] + p4 * x[6] * cos(x[6] - x[4]) / x[5] - p4 * atan(lv * x[3] / x[5]) * cos(x[6] - x[4]) / x[5] + p5 * x[6] * u[1] / x[5]
			- p5 * x[4] * u[1] / x[5] + p6 * atan(lh * x[3] / x[5]) / x[5] - p7 * x[4] * u[1] / x[5] + (p9 * x[3]) * (p9 * x[3]) * x[4] * x[5]) * delta_t
			+ x[4];
	float dx5 = (p4 * x[6] * x[4] - (p4 * x[6]) * (p4 * x[6]) - p4 * atan(lv * x[3] / x[5]) * x[6] + p4 * atan(lv * x[3] / x[5]) * x[4]
			+ p5 * u[1] * cos(x[6] - x[4]) + p6 * atan(lh * x[3] / x[5]) * x[4] + p7 * u[1] - p8 * x[5] - (p9 * x[5]) * (p9 * x[5]) * x[3] * x[3]
			- p10 * x[5] * x[5]) * delta_t + x[5];
	float dx6 = (-x[6] / p11 + u[0] / p11) * delta_t + x[6];



	x[0] = dx0;
	x[1] = dx1;
	x[2] = dx2;
	x[3] = dx3;
	x[4] = dx4;
	x[5] = dx5;
	x[6] = dx6;
	return x;
}

//float * f_cont(float *x, float t, float *u) {
//
////// Meine Modellgleichungen aus Matlab
//	float dx0 = x[5] * cos(x[2] + x[4]);
//	float dx1 = x[5] * sin(x[2] + x[4]);
//	float dx2 = x[3];
//	float dx3 = (-p1 * atan(lh * x[3] / x[5]) + p2 * x[6] * cos(x[6]) - p2 * atan(lv * x[3] / x[5]) * cos(x[6]) + p3 * u[1] * x[6]);
//	float dx4 = (-x[3] + p4 * x[6] * cos(x[6] - x[4]) / x[5] - p4 * atan(lv * x[3] / x[5]) * cos(x[6] - x[4]) / x[5] + p5 * x[6] * u[1] / x[5]
//			- p5 * x[4] * u[1] / x[5] + p6 * atan(lh * x[3] / x[5]) / x[5] - p7 * x[4] * u[1] / x[5] + square(p9 * x[3]) * x[4] * x[5]);
//	float dx5 = (p4 * x[6] * x[4] - square(p4 * x[6]) - p4 * atan(lv * x[3] / x[5]) * x[6] + p4 * atan(lv * x[3] / x[5]) * x[4] + p5 * u[1] * cos(x[6] - x[4])
//			+ p6 * atan(lh * x[3] / x[5]) * x[4] + p7 * u[1] - p8 * x[5] - square(p9 * x[5]) * square(x[3]) - p10 * square(x[5]));
//// dx5 = 0
//	float ret[] = { dx0, dx1, dx2, dx3, dx4, dx5 };
//	return ret;
//}
//
//float * f_cont_t_first(float t, float *x, float *u) {
//
////// Meine Modellgleichungen aus Matlab
//	float dx0 = x[5] * cos(x[2] + x[4]);
//	float dx1 = x[5] * sin(x[2] + x[4]);
//	float dx2 = x[3];
//	float dx3 = (-p1 * atan(lh * x[3] / x[5]) + p2 * x[6] * cos(x[6]) - p2 * atan(lv * x[3] / x[5]) * cos(x[6]) + p3 * u[1] * x[6]);
//	float dx4 = (-x[3] + p4 * x[6] * cos(x[6] - x[4]) / x[5] - p4 * atan(lv * x[3] / x[5]) * cos(x[6] - x[4]) / x[5] + p5 * x[6] * u[1] / x[5]
//			- p5 * x[4] * u[1] / x[5] + p6 * atan(lh * x[3] / x[5]) / x[5] - p7 * x[4] * u[1] / x[5] + square(p9 * x[3]) * x[4] * x[5]);
//	float dx5 = (p4 * x[6] * x[4] - square(p4 * x[6]) - p4 * atan(lv * x[3] / x[5]) * x[6] + p4 * atan(lv * x[3] / x[5]) * x[4] + p5 * u[1] * cos(x[6] - x[4])
//			+ p6 * atan(lh * x[3] / x[5]) * x[4] + p7 * u[1] - p8 * x[5] - square(p9 * x[5]) * x[3] * x[3] - p10 * x[5] * x[5]);
//// dx5 = 0
//	float dx6 = -x[6] / p11 + u[0] / p11;
//	float ret[] = { dx0, dx1, dx2, dx3, dx4, dx5, dx6 };
//	return ret;
//}

