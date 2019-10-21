/**
 * @file   	%<%NAME%>%.%<%EXTENSION%>%
 * @author 	%<%USER%>%
 * @version
 * @date 	%<%DATE%>%, %<%TIME%>%
 * @brief
 */


/*#### |Begin| --> Секция - "Include" ########################################*/
#include "Lib_A_PGCS_pos_gnss_correct_system.h"
/*#### |End  | <-- Секция - "Include" ########################################*/


/*#### |Begin| --> Секция - "Глобальные переменные" ##########################*/
/*#### |End  | <-- Секция - "Глобальные переменные" ##########################*/


/*#### |Begin| --> Секция - "Локальные переменные" ###########################*/
/*#### |End  | <-- Секция - "Локальные переменные" ###########################*/


/*#### |Begin| --> Секция - "Прототипы локальных функций" ####################*/
/*#### |End  | <-- Секция - "Прототипы локальных функций" ####################*/


/*#### |Begin| --> Секция - "Описание глобальных функций" ####################*/

/* @todo Эти определения перенеси в .h и сделай так, чтобы компилятор "ругался", 
 * если явно не выбран метод проекции */
#ifndef __PGCS_BACKPROJECTMETHOD
	#define __PGCS_BackProjectCoordSys(x)	PGCS_FlatToLLA(x)
#elif (__PGCS_BACKPROJECTMETHOD == 1)
	#define __PGCS_BackProjectCoordSys(x)	PGCS_FlatToLLA(x)
#elif (__PGCS_BACKPROJECTMETHOD == 2)
	#define __PGCS_BackProjectCoordSys(x)	PGCS_ECEFToLLAAdd(x)
#endif

void
PGCS_StructInit(
	pgcs_data_init_s *pInit_s)
{
	/* Сброс скалярных параметров в значения по умолчанию */
	pInit_s->scalParams_s.alpha 	= (__PGCS_FPT__) 1.0;
	pInit_s->scalParams_s.beta 		= (__PGCS_FPT__) 2.0;
	pInit_s->scalParams_s.kappa 	= (__PGCS_FPT__) 0.0;

	/* Сброс периода интегрирования */
	pInit_s->dt = (__PGCS_FPT__) 0.0;

	/* Сброс в нуль: */
	for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)
	{
		pInit_s->Q_mat_a[i] = (__PGCS_FPT__) 0.0;
		pInit_s->R_mat_a[i] = (__PGCS_FPT__) 0.0;
		pInit_s->state_a[i] = (__PGCS_FPT__) 0.0;
	}
}

void
PGCS_Init_All(
	pgcs_data_s *pData_s,
	pgcs_data_init_s *pInit_s)
{
	/* @todo Сделанй проверку периода интегрирования, если он равен нулю, чтобы эта функиця зацикливалась */

	/* Обновление периода интегрирования */
	PGCS_UpdateDt(pData_s, pInit_s->dt);

	/* Задание указателей вектора dpos на результат интегрирования в структурах типа ninteg_trapz_s */
	for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)

		/* Есть предположение, что не до конца ясен способ работы с библиотекой интегрирования */
		pData_s->kinData_s.flat_dpos[i] = &(pData_s->kinData_s.flat_pos_integ[i].deltaData);

	/* @todo А где ты инициализируешь структуру типа "ninteg_trapz_s", иначе функция интегрирования работать не будет */
}

void
PGCS_UpdatePosState(
	pgcs_data_s *pData_s)
{
	if __PGCS_IsFlagVelDataUpdateSet()
	{
		PGCS_IntegrateLLA(pData_s);
		__PGCS_ReSetFlagVelDataUpdate();
	}

	/* Выполнение проекции приращения местоположения из нормальной Земной СК 
	 * (модели плоской Земли) в приращение долготы/широты/высоты */
	__PGCS_BackProjectCoordSys(pData_s);

}

void
PGCS_CopyVelInWorldFrame(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pVel)
{
	if !(__PGCS_IsFlagVelDataUpdateSet())
	{
		/* @todo раскрути цикл, количество параметров не измениться */
		for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)
			pData_s->kinData_s.flat_vel[i] = *pVel++;
		__PGCS_SetFlagVelDataUpdate();
	}
}

void
PGCS_CopyLatLonAltInWorldFrame(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pLatLonAlt)
{
	if !(__PGCS_IsFlagPosDataUpdateSet())
	{
		/* @todo раскрути цикл, количество параметров не измениться */
		for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)
			pData_s->kinData_s.lla_pos_zero[i] = *pLatLonAlt++;
		__PGCS_SetFlagPosDataUpdate();
	}
}
/*#### |End  | <-- Секция - "Описание глобальных функций" ####################*/


/*#### |Begin| --> Секция - "Описание локальных функций" #####################*/

/*void
PGCS_ECEFToLLA(
	pgcs_data_s *pData_s)
{
	//ToDo

	__PGCS_FPT__ a = 6378137.0;            earth semimajor axis in meters
	__PGCS_FPT__ f = 1. / 298.257223563;   reciprocal flattening
	__PGCS_FPT__ b = a * (1. - f);                semi-minor axis
	__PGCS_FPT__ b2 = b * b;

	__PGCS_FPT__ e2 = 2.*f - (f * f);             first eccentricity squared
	__PGCS_FPT__ ep2 = f * (2. - f) / ((1. - f) * (1. - f));  second eccentricity squared
	__PGCS_FPT__ E2 = a * a - b2;


	__PGCS_FPT__ z2 = pData_s->kinData_s.flat_pos_cur[2] * pData_s->kinData_s.flat_pos_cur[2];
	__PGCS_FPT__ r2 = pData_s->kinData_s.flat_pos_cur[0] * pData_s->kinData_s.flat_pos_cur[0] + pData_s->kinData_s.flat_pos_cur[1] * pData_s->kinData_s.flat_pos_cur[1];
	__PGCS_FPT__ r = __PGCS_sqrt(r2);
	__PGCS_FPT__ F = 54.*b2 * z2;
	__PGCS_FPT__ G = r2 + (1 - e2) * z2 - e2 * E2;
	__PGCS_FPT__ c = (e2 * e2 * F * r2) / (G * G * G);
	__PGCS_FPT__ s = __PGCS_pow((1 + c + __PGCS_sqrt(c * c + 2 * c)), 1. / 3.);
	__PGCS_FPT__ s1 = 1 + s + 1 / s;
	__PGCS_FPT__ P = F / (3 * s1 * s1 * G * G);
	__PGCS_FPT__ Q = __PGCS_sqrt(1 + 2 * e2 * e2 * P);
	__PGCS_FPT__ ro = -(e2 * P * r) / (1 + Q) + __PGCS_sqrt((a * a / 2) * (1 + 1 / Q) - ((1 - e2) * P * z2) / (Q *
                   (1 + Q)) - P * r2 / 2);
	__PGCS_FPT__ tmp = (r - e2 * ro) * (r - e2 * ro);
	__PGCS_FPT__ U = __PGCS_sqrt(tmp + z2);
	__PGCS_FPT__ V = __PGCS_sqrt(tmp + (1 - e2) * z2);
	__PGCS_FPT__ zo = (b2 * pData_s->kinData_s.flat_pos_cur[2]) / (a * V);

	pData_s->kinData_s.lla_pos_cur[0] = __PGCS_atan((pData_s->kinData_s.flat_pos_cur[2] + ep2 * zo) / r);
	pData_s->kinData_s.lla_pos_cur[1] = __PGCS_atan2(pData_s->kinData_s.flat_pos_cur[1], pData_s->kinData_s.flat_pos_cur[0]);
	pData_s->kinData_s.lla_pos_cur[2] = U * (1 - b2 / (a * V));



} */


/*void
PGCS_ECEFToLLAAdd(
	pgcs_data_s *pData_s)
{
	//ToDo

	PGCS_ECEFToLLA(pData_s);
	for (uint8_t i=0; i<3; i++)
		pData_s->kinData_s.lla_pos_cur[i] += pData_s->kinData_s.lla_pos_zero[i];
}*/

void
PGCS_FlatToLLA(
	pgcs_data_s *pData_s)
{
	/* @todo константу объяви как #define */
	__PGCS_FPT__ re = (__PGCS_FPT__) 6378137.0;
	__PGCS_FPT__ re_c = re * __PGCS_cos((pi / ((__PGCS_FPT__) 180.0) * __PGCS_fabs(pData_s->kinData_s.lla_pos_zero[0]));

	/* @todo Ты ведь выполняешь проекцию не местоположения, а приращения местоположения, может стоит это явно указать в имени переменной 
	 * pData_s->kinData_s.lla_pos[0] заменить на pData_s->kinData_s.lla_deltaPos[0]
	 */
	pData_s->kinData_s.lla_pos[0] = *pData_s->kinData_s.flat_dpos[1] * ((__PGCS_FPT__) 180.0) / (pi * re) + pData_s->kinData_s.lla_pos_zero[0];
	pData_s->kinData_s.lla_pos[1] = *pData_s->kinData_s.flat_dpos[0] * ((__PGCS_FPT__) 180.0) / (pi * re_c) + pData_s->kinData_s.lla_pos_zero[1];
	pData_s->kinData_s.lla_pos[2] = *pData_s->kinData_s.flat_dpos[2] + pData_s->kinData_s.lla_pos_zero[2];


	/* @todo смотри описание этого флага в .h */
	__PGCS_SetFlagPosDataUpdate();
}

/* @todo Эта функция ведь получает приращение местоположения в нормальной земной СК (т.е. в модели плоской Земли), почему функция
 * называется PGCS_IntegrateLLA() вместо PGCS_IntegrateFLatEath()?  */
void 
PGCS_IntegrateLLA(
	pgcs_data_s *pData_s)
{
	/* @todo раскрути цикл, количество параметров не измениться */
	for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)
	{
		/* Получение приращения местоположения в нормальной земной СК */
		NINTEG_Trapz(&(pData_s->kinData_s.flat_pos_integ[i]), pData_s->kinData_s.flat_vel[i]);
		//pData_s->kinData_s.flat_dpos[i] = pData_s->kinData_s.flat_pos_integ[i].deltaData;
	}

	__PGCS_SetFlagDPosDataUpdate();
}

void
PGCS_UpdateDt(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ dt)
{
	for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)
		pData_s->kinData_s.flat_pos_integ[i].dT = dt;
}

/*#### |End  | <-- Секция - "Описание локальных функций" #####################*/


/*#### |Begin| --> Секция - "Обработчики прерываний" #########################*/
/*#### |End  | <-- Секция - "Обработчики прерываний" #########################*/

/*############################################################################*/
/*############################ END OF FILE  ##################################*/
/*############################################################################*/
