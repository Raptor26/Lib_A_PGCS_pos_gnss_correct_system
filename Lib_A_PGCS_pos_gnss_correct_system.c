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
void
PGCS_ECEFToLLAAdd2(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pDPos);

void
PGCS_FlatToLLA2(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pDPos);

void
PGCS_FlatToLLA1(
	pgcs_data_s *pData_s);

void 
PGCS_IntegrateFlat(
	pgcs_data_s *pData_s);

void
PGCS_UpdateDt(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ dt);

void __PGCS_FNC_ONCE_MEMORY_LOCATION
PGSS_Init_MatrixStructs(
	pgcs_data_s 		*pData_s,
	ukfsif_all_data_s 	*pMatrixPointers_s);

/*#### |End  | <-- Секция - "Прототипы локальных функций" ####################*/


/*#### |Begin| --> Секция - "Описание глобальных функций" ####################*/

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

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      03-сен-2019
 *
 * @brief     Функция выполняет инициализацию всех параметров, необходимых
 *            для работы UKF, операций интегрирования
 *
 * @param[in, out] 	*pData_s:	Указатель на структуру данных со всеми параметрами
 * 
 * @param[in]   	*pInit_s:   Указатель на структуру с параметрами инициализации
 *
 */
void
PGCS_Init_All(
	pgcs_data_s *pData_s,
	pgcs_data_init_s *pInit_s)
{
	/* Инициализация интегральных структур */
	ninteg_trapz_init_s trapzInit_s;
	NINTEG_Trapz_StructInit(&trapzInit_s);
	trapzInit_s.accumulate_flag = NINTEG_DISABLE;
	trapzInit_s.integratePeriod = pInit_s->dt;

	for (uint8_t i = 0; i < PGCS_LEN_STATE; i++)
	{
		ninteg_fnc_status_e trapzInitStatus_e =
			NINTEG_Trapz_Init(
				&pData_s->kinData_s.flat_pos_integ[i],
				&trapzInit_s);
		/* Зависнуть если ошибка инициализации */
		while (trapzInitStatus_e != NINTEG_SUCCESS);
	}

	/* Инициализация всех структур матриц */
	PGSS_Init_MatrixStructs(
		pData_s,
		&pData_s->ukfData_s.ukfsifMatrixPointers_s);

	/* Вычисление корня квадратного из (lambda + len) и запись в поле структуры */
	pData_s->ukfData_s.scalar_s.sqrtLamLen =
		__PGCS_sqrt(
			UKFSIF_GetLambda(
				PGCS_LEN_STATE,
				pInit_s->scalParams_s.alpha,
				pInit_s->scalParams_s.kappa) + PGCS_LEN_STATE);

	/* Установка периода интегрирования */
	#if defined (__UKFMO_CHEKING_ENABLE__)
	if (pInit_s->dt == (__PGCS_FPT__)0.0)
	{
		__UKFMO_ALL_INTERRUPTS_DIS();
		while (1);
	}
	#endif

	/* Обновление периода интегрирования */
	__PGCS_UpdateDt(pData_s, pInit_s->dt);

	/* Инициализация вектора muMean */
	UKFSIF_InitWeightVectorMean(
		&pInit_s->scalParams_s,
		pData_s->ukfData_s.muMean_s.memForMatrix[0u],
		PGCS_LEN_STATE);

	/* Инициализация вектора muCov */
	UKFSIF_InitWeightVectorCov(
		&pInit_s->scalParams_s,
		pData_s->ukfData_s.muCovar_s.memForMatrix[0u],
		PGCS_LEN_STATE);

	/* Заполнение матрицы Q */
	PGCS_Init_NoiseMatrix(
		&pData_s->ukfData_s.noiseMatrix_s.QMat_s.mat_s,
		pInit_s->Q_mat_a);

	/* Заполнение матрицы R */
	PGCS_Init_NoiseMatrix(
		&pData_s->ukfData_s.noiseMatrix_s.RMat_s.mat_s,
		pInit_s->R_mat_a);

	/* Заполнение матрицы P */
	UKFMO_MatrixIdentity(
		&pData_s->ukfData_s.P_predict_s.mat_s);

	/* Инициализация вектора пространства состояний */
	for (size_t i = 0u; i < __UKFMO_GetRowNumb(&pData_s->ukfData_s.x_posteriori_s.mat_s); i++)
	{
		/* Получить индекс ячейки массива */
		size_t idx =
			__UKFMO_GetIndexInOneFromTwoDim(&pData_s->ukfData_s.x_posteriori_s.mat_s, i, 0u);
		pData_s->ukfData_s.x_posteriori_s.mat_s.pData[idx]
			= pInit_s->state_a[i];
	}

}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция производит обновление позиции путём выполнения операций
 * 			  интегрирования, а также всех этапов работы UKF
 *
 * @param[in, out]    *pData_s:    Указатель на структуру данных, содержащую
 * 								   данные кинематики и UKF
 */
void
PGCS_UpdatePosState(
	pgcs_data_s *pData_s)
{
	if __PGCS_IsFlagVelDataUpdateSet()
	{
		PGCS_IntegrateFlat(pData_s);
		__PGCS_ReSetFlagVelDataUpdate();
	}

	/* Выполнение проекции приращения местоположения из нормальной Земной СК 
	 * (модели плоской Земли) в приращение долготы/широты/высоты */
	__PGCS_BackProjectCoordSys1(pData_s);

}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция производит запись вектора скорости во внутреннюю структуру
 *
 * @param[in, out]    *pData_s:    Указатель на структуру данных, содержащую
 * 								   данные кинематики
 * 								   
 * @param[in]    	  *pVel:       Указатель на массив (вектор) скорости по трём осям
 */
void
PGCS_SetCurrentFlatVelocity(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pVel)
{
	if !(__PGCS_IsFlagVelDataUpdateSet())
	{
		pData_s->kinData_s.flat_vel[0] = *pVel++;
		pData_s->kinData_s.flat_vel[1] = *pVel++;
		pData_s->kinData_s.flat_vel[2] = *pVel;
		__PGCS_SetFlagVelDataUpdate();
	}
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция производит запись координаты, полученной с GNSS,
 * 		      во внутреннюю структуру
 *
 * @param[in, out]    *pData_s:    		Указатель на структуру данных, содержащую
 * 								   		данные кинематики
 * 								   
 * @param[in]    	  *pLatLonAlt:      Указатель на массив (вектор) координаты ДШВ
 */
void
PGCS_SetCurrentLLAPos(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pLatLonAlt)
{
	if !(__PGCS_IsFlagPosDataUpdateSet())
	{
		pData_s->kinData_s.lla_pos_gnss[0] = *pLatLonAlt++;
		pData_s->kinData_s.lla_pos_gnss[1] = *pLatLonAlt++;
		pData_s->kinData_s.lla_pos_gnss[2] = *pLatLonAlt;
		__PGCS_SetFlagPosDataUpdate();
	}
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция производит запись опорной нулевой координаты, полученной
 * 			  с GNSS, во внутреннюю структуру. Используется для корректной работы
 * 			  функций обратной проекции.
 *
 * @param[in, out]    *pData_s:    		Указатель на структуру данных, содержащую
 * 								   		данные кинематики
 * 								   
 * @param[in]    	  *pLatLonAlt:      Указатель на массив (вектор) нулевой
 * 										координаты ДШВ
 */
void
PGCS_SetZeroLLAPos(
  pgcs_data_s *pData_s,
  __PGCS_FPT__ *pLatLonAlt)
{
		pData_s->kinData_s.lla_pos_zero[0] = *pLatLonAlt++;
		pData_s->kinData_s.lla_pos_zero[1] = *pLatLonAlt++;
		pData_s->kinData_s.lla_pos_zero[2] = *pLatLonAlt;
}
/*#### |End  | <-- Секция - "Описание глобальных функций" ####################*/


/*#### |Begin| --> Секция - "Описание локальных функций" #####################*/

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция обратного проецирования координат из прямоугольной геоцентрической
 * 			  системы (ECEF) во всемирную систему (WGS84)
 *
 * @param[in,out]	*pData_s:    Указатель на структуру данных, в которой содержатся
 * 								 кинетические данные
 * 								 
 * @param[in]    	*pDPos:      Указатель на массив (вектор) координат приращения позиции
 * 								 в проекционной системе
 */
void
PGCS_ECEFToLLAAdd2(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pDPos)
{
	/*					ToDo						*/
	__PGCS_FPT__ f = 1. / 298.257223563;		eciprocal flattening
	__PGCS_FPT__ b = PGCS_RE * (1. - f);		semi-minor axis
	__PGCS_FPT__ b2 = b * b;

	__PGCS_FPT__ e2 = 2.*f - (f * f);							first eccentricity squared
	__PGCS_FPT__ ep2 = f * (2. - f) / ((1. - f) * (1. - f));	second eccentricity squared
	__PGCS_FPT__ E2 = PGCS_RE * PGCS_RE - b2;


	__PGCS_FPT__ z2 = pData_s->kinData_s.flat_dpos[2] * pData_s->kinData_s.flat_dpos[2];
	__PGCS_FPT__ r2 = pData_s->kinData_s.flat_dpos[0] * pData_s->kinData_s.flat_dpos[0] + pData_s->kinData_s.flat_dpos[1] * pData_s->kinData_s.flat_dpos[1];
	__PGCS_FPT__ r = __PGCS_sqrt(r2);
	__PGCS_FPT__ F = 54.*b2 * z2;
	__PGCS_FPT__ G = r2 + (1 - e2) * z2 - e2 * E2;
	__PGCS_FPT__ c = (e2 * e2 * F * r2) / (G * G * G);
	__PGCS_FPT__ s = __PGCS_pow((1 + c + __PGCS_sqrt(c * c + 2 * c)), 1. / 3.);
	__PGCS_FPT__ s1 = 1 + s + 1 / s;
	__PGCS_FPT__ P = F / (3 * s1 * s1 * G * G);
	__PGCS_FPT__ Q = __PGCS_sqrt(1 + 2 * e2 * e2 * P);
	__PGCS_FPT__ ro = -(e2 * P * r) / (1 + Q) + __PGCS_sqrt((PGCS_RE * PGCS_RE / 2) * (1 + 1 / Q) - ((1 - e2) * P * z2) / (Q *
                   (1 + Q)) - P * r2 / 2);
	__PGCS_FPT__ tmp = (r - e2 * ro) * (r - e2 * ro);
	__PGCS_FPT__ U = __PGCS_sqrt(tmp + z2);
	__PGCS_FPT__ V = __PGCS_sqrt(tmp + (1 - e2) * z2);
	__PGCS_FPT__ zo = (b2 * pData_s->kinData_s.flat_dpos[2]) / (PGCS_RE * V);

	pData_s->kinData_s.lla_pos[0] = __PGCS_atan((pData_s->kinData_s.flat_dpos[2] + ep2 * zo) / r) + pData_s->kinData_s.lla_pos_zero[0];
	pData_s->kinData_s.lla_pos[1] = __PGCS_atan2(pData_s->kinData_s.flat_dpos[1], pData_s->kinData_s.flat_dpos[0]) + pData_s->kinData_s.lla_pos_zero[1];
	pData_s->kinData_s.lla_pos[2] = U * (1 - b2 / (PGCS_RE * V)) + pData_s->kinData_s.lla_pos_zero[2];
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция обратного проецирования координат из прямоугольной геоцентрической
 * 			  системы (ECEF) во всемирную систему (WGS84)
 *
 * @param[in,out]   	*pData_s:    Указатель на структуру данных, в которой содержатся
 * 								 	 кинетические данные
 */
void
PGCS_ECEFToLLAAdd1(
	pgcs_data_s *pData_s)
{
	PGCS_ECEFToLLAAdd2(pData_s, pData_s->kinData_s.flat_dpos);
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция обратного проецирования координат из плоскоземельной
 * 			  системы (ECEF) во всемирную систему (WGS84)
 *
 * @param[in,out]		*pData_s:    Указатель на структуру данных, в которой содержатся
 * 								 	 кинетические данные
 * 								 
 * @param[in]			*pDPos:      Указатель на массив (вектор) координат приращения позиции
 * 								 	 в проекционной системе
 */
void
PGCS_FlatToLLA2(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ *pDPos)
{
	__PGCS_FPT__ re_c = PGCS_RE * __PGCS_cos((PGCS_PI / ((__PGCS_FPT__) 180.0) * __PGCS_fabs(pData_s->kinData_s.lla_pos_zero[0]));

	pData_s->kinData_s.lla_pos[0] = *(pDPos + 1) * ((__PGCS_FPT__) 180.0) / (PGCS_PI * PGCS_RE) + pData_s->kinData_s.lla_pos_zero[0];
	pData_s->kinData_s.lla_pos[1] = *(pDPos + 0) * ((__PGCS_FPT__) 180.0) / (PGCS_PI * re_c) + pData_s->kinData_s.lla_pos_zero[1];
	pData_s->kinData_s.lla_pos[2] = *(pDPos + 2) + pData_s->kinData_s.lla_pos_zero[2];
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция обратного проецирования координат из плоскоземельной
 * 			  системы (ECEF) во всемирную систему (WGS84)
 *
 * @param[in,out]   	*pData_s:    Указатель на структуру данных, в которой содержатся
 * 								     кинетические данные
 */
void
PGCS_FlatToLLA1(
	pgcs_data_s *pData_s)
{
	PGCS_FlatToLLA2(pData_s, pData_s->kinData_s.flat_dpos);
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief     Функция, производящая интегрирование вектора скоростей и заполняющая
 * 			  вектор приращения координат в проекционной системе последним хранимым
 * 			  значением
 * 			  
 *
 * @param[in,out]		*pData_s:    Указатель на структуру данных, в которой содержатся
 * 						    	 	 кинетические данные
 */
void 
PGCS_IntegrateFlat(
	pgcs_data_s *pData_s)
{
		NINTEG_Trapz(&(pData_s->kinData_s.flat_pos_integ[0]), pData_s->kinData_s.flat_vel[0]);
		NINTEG_Trapz(&(pData_s->kinData_s.flat_pos_integ[1]), pData_s->kinData_s.flat_vel[1]);
		NINTEG_Trapz(&(pData_s->kinData_s.flat_pos_integ[2]), pData_s->kinData_s.flat_vel[2]);

		pData_s->kinData_s.flat_dpos[0] = NINTEG_TrapzGetLastVal(pData_s->kinData_s.flat_pos_integ[0]);
		pData_s->kinData_s.flat_dpos[1] = NINTEG_TrapzGetLastVal(pData_s->kinData_s.flat_pos_integ[1]);
		pData_s->kinData_s.flat_dpos[2] = NINTEG_TrapzGetLastVal(pData_s->kinData_s.flat_pos_integ[2]);
}

/*-------------------------------------------------------------------------*//**
 * @author    Konstantin Ganshin
 * @date      29-Oct-2019
 *
 * @brief    Функция производит обновление периода интегрирования в структурах,
 * 			 используемых для операций интегрирования
 *
 * @param[in,out]		*pData_s:    Указатель на структуру данных, в которой содержатся
 * 						   	  	     кинетические данные
 * 						   	  	   
 * @param[in]				  dt:    Новый период интегрирования
 */
void
PGCS_UpdateDt(
	pgcs_data_s *pData_s,
	__PGCS_FPT__ dt)
{
	pData_s->kinData_s.flat_pos_integ[0].dT = dt;
	pData_s->kinData_s.flat_pos_integ[1].dT = dt;
	pData_s->kinData_s.flat_pos_integ[2].dT = dt;
}

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @author    Konstantin Ganshin
 * @date      29-сен-2019
 *
 * @brief    Инициализация указателей на области памяти, в которых
 *           содержатся матрицы
 *
 * @param[out] 	*pData_s:			Указатель на структуру данных, в которой содержаться
 * 									параметры, необходимые для работы UKF
 * 									
 * @param[in]  	*pMatrixPointers_s: Указатель на структуру данных, содержащую 
 * 									указатели на области памяти матричных структур
 */
void __PGCS_FNC_ONCE_MEMORY_LOCATION
PGSS_Init_MatrixStructs(
	pgcs_data_s 		*pData_s,
	ukfsif_all_data_s 	*pMatrixPointers_s)
{
	/* Объявление структуры для инициализации указателей на матрицы */
	ukfsif_all_data_init_s    initMatrixPointers_s;
	UKFIS_StructInit		(&initMatrixPointers_s);

	/* Инициализация матрицы шума Q */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.noiseMatrix_s.QMat_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.ukfData_s.noiseMatrix_s.QMat_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.noiseMatrix_s.QMat_s.mat_s,
		sizeof(pData_s->ukfData_s.noiseMatrix_s.QMat_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_Q] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.noiseMatrix_s.QMat_s.mat_s);

	/* Инициализация матрицы шума R */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.noiseMatrix_s.RMat_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.noiseMatrix_s.RMat_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.noiseMatrix_s.RMat_s.mat_s,
		sizeof(pData_s->ukfData_s.noiseMatrix_s.RMat_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_R] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.noiseMatrix_s.RMat_s.mat_s);

	/* Инициализация матрицы пространства состояний */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.stateMat_s.mat_s, 			/* !< Указатель на структуру матрицы */
		PGCS_LEN_MATRIX_ROW, 					/* !< Количество строк */
		PGCS_LEN_MATRIX_COL,					/* !< Количество столбцов */
		pData_s->ukfData_s.stateMat_s.memForMatrix[0u] 	/* !< Указатель на область памяти для хранения матрицы */
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.stateMat_s.mat_s,
		sizeof (pData_s->ukfData_s.stateMat_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_x_LxL] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.stateMat_s.mat_s);

	/* Инициализация вектора пространства состояний */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.x_apriori_s.mat_s,
		PGCS_LEN_STATE,
		1u,
		pData_s->ukfData_s.x_apriori_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.x_apriori_s.mat_s,
		sizeof(pData_s->ukfData_s.x_apriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_x_apriori] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.x_apriori_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.x_posteriori_s.mat_s,
		PGCS_LEN_STATE,
		1u,
		pData_s->ukfData_s.x_posteriori_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.x_posteriori_s.mat_s,
		sizeof(pData_s->ukfData_s.x_posteriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_x_posteriori] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.x_posteriori_s.mat_s);

	/* Инициализация матрицы сигма-точек */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.chiSigmaMat_s.mat_s,
		PGCS_LEN_SIGMA_ROW,
		PGCS_LEN_SIGMA_COL,
		pData_s->ukfData_s.chiSigmaMat_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.chiSigmaMat_s.mat_s,
		sizeof(pData_s->ukfData_s.chiSigmaMat_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_chi_predict] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.chiSigmaMat_s.mat_s);

	/* Инициализация матрицы сигма-точек (после функции преобразования) */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.chiSigmaPostMat_s.mat_s,
		PGCS_LEN_SIGMA_ROW,
		PGCS_LEN_SIGMA_COL,
		pData_s->ukfData_s.chiSigmaPostMat_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)	&pData_s->ukfData_s.chiSigmaPostMat_s.mat_s,
		sizeof				(pData_s->ukfData_s.chiSigmaPostMat_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_chi_apriori] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.chiSigmaPostMat_s.mat_s);

	/* Инициализация матрицы квадратного корня от матрицы ковариации "P" */
	UKFMO_MatrixInit(
		&pData_s->ukfData_s.sqrtP_apriori_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.sqrtP_apriori_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.sqrtP_apriori_s.mat_s,
		sizeof(pData_s->ukfData_s.sqrtP_apriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_P_sqrt] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.sqrtP_apriori_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.muMean_s.mat_s,
		PGCS_LEN_SIGMA_COL,
		1u,
		pData_s->ukfData_s.muMean_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.muMean_s.mat_s,
		sizeof(pData_s->ukfData_s.muMean_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_muMean] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.muMean_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.muCovar_s.mat_s,
		PGCS_LEN_SIGMA_COL,
		1u,
		pData_s->ukfData_s.muCovar_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.muCovar_s.mat_s,
		sizeof(pData_s->ukfData_s.muCovar_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_muCovar] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.muCovar_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.chi_apriory_minus_x_apriory_s.mat_s,
		PGCS_LEN_STATE,
		1u,
		pData_s->ukfData_s.chi_apriory_minus_x_apriory_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*)&pData_s->ukfData_s.chi_apriory_minus_x_apriory_s.mat_s,
		sizeof (pData_s->ukfData_s.chi_apriory_minus_x_apriory_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.chi_apriory_minus_x_apriory_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.chi_apriory_minus_x_apriory_Transpose_s.mat_s,
		1u,
		PGCS_LEN_STATE,
		pData_s->ukfData_s.chi_apriory_minus_x_apriory_Transpose_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		(ukfmo_matrix_s*) &pData_s->ukfData_s.chi_apriory_minus_x_apriory_Transpose_s.mat_s,
		sizeof(pData_s->ukfData_s.chi_apriory_minus_x_apriory_Transpose_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_chi_priory_MINUS_x_priory_TRANSPOSE] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.chi_apriory_minus_x_apriory_Transpose_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.resultOfMult2Matrix_s.mat_s,
		PGCS_LEN_STATE,
		PGCS_LEN_STATE,
		pData_s->ukfData_s.resultOfMult2Matrix_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.resultOfMult2Matrix_s.mat_s,
		sizeof(pData_s->ukfData_s.resultOfMult2Matrix_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_result_of_mult_2_matrix] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.resultOfMult2Matrix_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.P_apriori_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.P_apriori_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.P_apriori_s.mat_s,
		sizeof(pData_s->ukfData_s.P_apriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_P_apriory] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.P_apriori_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.psi_apriori_s.mat_s,
		PGCS_LEN_SIGMA_ROW,
		PGCS_LEN_SIGMA_COL,
		pData_s->ukfData_s.psi_apriori_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.psi_apriori_s.mat_s,
		sizeof(pData_s->ukfData_s.psi_apriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_psi_apriori] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.psi_apriori_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.y_apriori_s.mat_s,
		PGCS_LEN_STATE,
		1u,
		pData_s->ukfData_s.y_apriori_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.y_apriori_s.mat_s,
		sizeof(pData_s->ukfData_s.y_apriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_y_apriori] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.y_apriori_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.Pyy_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.Pyy_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.Pyy_s.mat_s,
		sizeof(pData_s->ukfData_s.Pyy_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_Pyy] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.Pyy_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.psi_priory_MINUS_y_priory_TRANSPOSE.mat_s,
		1u,
		PGCS_LEN_STATE,
		pData_s->ukfData_s.psi_priory_MINUS_y_priory_TRANSPOSE.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.psi_priory_MINUS_y_priory_TRANSPOSE.mat_s,
		sizeof(pData_s->ukfData_s.psi_priory_MINUS_y_priory_TRANSPOSE.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_psi_priory_MINUS_y_priory_TRANSPOSE] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.psi_priory_MINUS_y_priory_TRANSPOSE.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.Pxy_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.Pxy_s.memForMatrix[0u]
	);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.Pxy_s.mat_s,
		sizeof(pData_s->ukfData_s.Pxy_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_Pxy] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.Pxy_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.PyyInv_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.PyyInv_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.PyyInv_s.mat_s,
		sizeof(pData_s->ukfData_s.PyyInv_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_Pyy_INV] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.PyyInv_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.K_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.K_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.K_s.mat_s,
		sizeof(pData_s->ukfData_s.K_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_K] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.K_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.y_posteriori_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		1u,
		pData_s->ukfData_s.y_posteriori_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.y_posteriori_s.mat_s,
		sizeof(pData_s->ukfData_s.y_posteriori_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_y_posteriori] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.y_posteriori_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.innovation_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		1u,
		pData_s->ukfData_s.innovation_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.innovation_s.mat_s,
		sizeof(pData_s->ukfData_s.innovation_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_innovation] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.innovation_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.P_predict_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.P_predict_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.P_predict_s.mat_s,
		sizeof(pData_s->ukfData_s.P_predict_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_P] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.P_predict_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.K_Transpose_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.K_Transpose_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.K_Transpose_s.mat_s,
		sizeof(pData_s->ukfData_s.K_Transpose_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_K_TRANSPOSE] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.K_Transpose_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.x_predict_temp_s.mat_s,
		PGCS_LEN_MATRIX_ROW,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.x_predict_temp_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.x_predict_temp_s.mat_s,
		sizeof(pData_s->x_predict_temp_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_x_LxL_TEMP] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.x_predict_temp_s.mat_s);

	UKFMO_MatrixInit(
		&pData_s->ukfData_s.x_predict_temp_ones_s.mat_s,
		1u,
		PGCS_LEN_MATRIX_COL,
		pData_s->ukfData_s.x_predict_temp_ones_s.memForMatrix[0u]);
	__UKFMO_CheckMatrixSize(
		&pData_s->ukfData_s.x_predict_temp_ones_s.mat_s,
		sizeof(pData_s->ukfData_s.x_predict_temp_ones_s.memForMatrix));
	initMatrixPointers_s.pMatrix_s_a[UKFSIF_INIT_x_1xL_ones_TEMP] =
		__PGCS_CheckMatrixStructValidation(
			&pData_s->ukfData_s.x_predict_temp_ones_s.mat_s);

	/* Копирование указателей на структуры матриц (Эта функция должна быть
	 * вызвана в конце) */
	UKFSIF_Init_SetMatrixPointers(
		pMatrixPointers_s,
		&initMatrixPointers_s,
		(uint16_t) PGCS_LEN_STATE);

}

static void __PGCS_FNC_ONCE_MEMORY_LOCATION
PGCS_Init_NoiseMatrix(
	ukfmo_matrix_s 	*pNoiseMat,
	__PGCS_FPT__ 	*pNoiseMatDiag)
{
	/* Сброс матрицы шумов в нуль */
	UKFMO_MatrixZeros(pNoiseMat);

	size_t i;
	for (i = 0u; i < pNoiseMat->numCols; i++)
	{
		if (*pNoiseMatDiag == ((__PGCS_FPT__) 0.0))
		{
			/* Если попоали сюда, значит диагональ матрицы шума не
			 * инициализирована */
			while (1);
		}
		pNoiseMat->pData[__UKFMO_GetIndexInOneFromTwoDim(pNoiseMat, i, i)] =
			*pNoiseMatDiag++;
	}
}

/*#### |End  | <-- Секция - "Описание локальных функций" #####################*/


/*#### |Begin| --> Секция - "Обработчики прерываний" #########################*/
/*#### |End  | <-- Секция - "Обработчики прерываний" #########################*/

/*############################################################################*/
/*############################ END OF FILE  ##################################*/
/*############################################################################*/
