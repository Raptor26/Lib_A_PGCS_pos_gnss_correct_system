/**
 * @file    %<%NAME%>%.%<%EXTENSION%>%
 * @author  %<%USER%>%
 * @version
 * @date    %<%DATE%>%, %<%TIME%>%
 * @brief
 */


#ifndef LIB_A_PGCS
#define LIB_A_PGCS


/*#### |Begin| --> Секция - "Include" ########################################*/
/*==== |Begin| --> Секция - "C libraries" ====================================*/
/*==== |End  | <-- Секция - "C libraries" ====================================*/

/*==== |Begin| --> Секция - "RTOS libraries ==================================*/
/*==== |End  | <-- Секция - "RTOS libraries ==================================*/

/*==== |Begin| --> Секция - "MK peripheral libraries" ========================*/
/*==== |End  | <-- Секция - "MK peripheral libraries" ========================*/

/*==== |Begin| --> Секция - "Extern libraries" ===============================*/
#include <math.h>
#include "Lib_A_UKFMO_ukf_matrix_operations.h"
#include "Lib_A_UKFSIF_ukf_standart_init_fnc.h"
#include "Lib_A_NINTEG_numerical_integration.h"
/*==== |End  | <-- Секция - "Extern libraries" ===============================*/
/*#### |End  | <-- Секция - "Include" ########################################*/


/*#### |Begin| --> Секция - "Определение констант" ###########################*/
#ifdef __PGCS_EXTERN_MODE_ENABLE__
	#include "macros_definitions.h"
#endif

/*==== |Begin| --> Секция определения метода обратной проекции ===============*/
#ifndef __PGCS_BACKPROJECTMETHOD
	#error "__PGCS_BACKPROJECTMETHOD isn't set. You must choose one of the backprojection methods and set it in macro list."
#elif (__PGCS_BACKPROJECTMETHOD == 1)
	#define __PGCS_BackProjectCoordSys2(x, y)	PGCS_FlatToLLA2(x, y)
	#define __PGCS_BackProjectCoordSys1(x)		PGCS_FlatToLLA1(x)
#elif (__PGCS_BACKPROJECTMETHOD == 2)
	#define __PGCS_BackProjectCoordSys2(x, y)	PGCS_ECEFToLLAAdd2(x, y)
	#define __PGCS_BackProjectCoordSys1(x)		PGCS_ECEFToLLAAdd1(x)
#else
	#error "__PGCS_BACKPROJECTMETHOD isn't correct. You must choose one of the backprojection methods and set it in macro list."
#endif

/*==== |End  | <-- Секция определения метода обратной проекции ===============*/

/*==== |Begin| --> Секция определения типа числа с плавающей точкой ==========*/
#ifndef __PGCS_FPT__
	#error "Please, set __PGCS_FPT__ to float or double in macros list"
#endif

#ifndef __PGCS_FPT_SIZE__
	#error "Please, set __PGCS_FPT_SIZE__ to 4 (that means float) or 8 (that means double) in macros list"
#endif

#if     __PGCS_FPT_SIZE__ == 4
	#define __PGCS_sqrt(x)		sqrtf(x)
	#define __PGCS_pow(x, y)	powf(x, y)
	#define __PGCS_atan(x)		atanf(x)
	#define __PGCS_atan2(x, y)	atan2f(x, y)
	#define __PGCS_cos(x)		cosf(x)
	#define __PGCS_fabs(x)		fabsf(x)

#elif   __PGCS_FPT_SIZE__ == 8
	#define __PGCS_sqrt(x)		sqrt(x)
	#define __PGCS_pow(x,y)		pow(x,y)
	#define __PGCS_atan(x)		atan(x)
	#define __PGCS_atan2(x, y)	atan2(x, y)
	#define __PGCS_cos(x)		cos(x)
	#define __PGCS_fabs(x)		fabs(x)
#else
	#error "Your compiler uses a non-standard floating point size"
#endif
/*==== |End  | <-- Секция определения типа числа с плавающей точкой ==========*/

/*==== |Begin| --> Секция - Макросы для встраиваемых функций =================*/
#ifdef __GNUC__

/* inline*/
#ifndef __PGCS_INLINE
	#define __PGCS_INLINE 			inline
#endif

/* static inline */
#ifndef __PGCS_STATIC_INLINE
	#define __PGCS_STATIC_INLINE	static inline
#endif

/* always inline */
#ifndef __PGCS_ALWAYS_INLINE
	#define __PGCS_ALWAYS_INLINE    inline __attribute__((always_inline)) static
#endif

/* force inline */
#ifndef __PGCS_FORCE_INLINE
	#define __PGCS_FORCE_INLINE 	inline __attribute__((always_inline))
#endif

#else
	#define __PGCS_INLINE
	#define __PGCS_STATIC_INLINE 		static
	#define __PGCS_ALWAYS_INLINE
#endif
/*==== |End  | <-- Секция - Макросы для встраиваемых функций =================*/

/*==== |Begin| --> Секция - Расположение функций библиотеки в специальной
 *                          области памяти ===================================*/
#ifdef __PGCS_FNC_ONCE_MEMORY_LOCATION_NAME__
	#ifdef __GNUC__
		#define __PGCS_FNC_ONCE_MEMORY_LOCATION  __attribute__ ((section(__PGCS_FNC_ONCE_MEMORY_LOCATION_NAME__)))
	#else
		#error "You defined the name of the memory area for the function location, but the type of your compiler is not supported by the library. You can delete the macro definition __PGCS_FNC_ONCE_MEMORY_LOCATION_NAME__ or extend the macro definition __PGCS_FNC_ONCE_MEMORY_LOCATION for your compiler type"
	#endif
#else
	#define __PGCS_FNC_ONCE_MEMORY_LOCATION
#endif

#ifdef __PGCS_FNC_LOOP_MEMORY_LOCATION_NAME__
	#ifdef (__GNUC__)
		#define __PGCS_FNC_LOOP_MEMORY_LOCATION  __attribute__ ((section(__PGCS_FNC_LOOP_MEMORY_LOCATION_NAME__)))
	#else
		#error "You defined the name of the memory area for the function location, but the type of your compiler is not supported by the library. You can delete the macro definition __PGCS_FNC_LOOP_MEMORY_LOCATION_NAME__ or extend the macro definition __PGCS_FNC_LOOP_MEMORY_LOCATION for your compiler type"
	#endif
#else
	#define __PGCS_FNC_LOOP_MEMORY_LOCATION
#endif
/*==== |End  | <-- Секция - Расположение функций библиотеки в специальной
 *                          области памяти ===================================*/

/*==== |Begin| --> Секция - Выбор метода численного интегрирования ===========*/
/*#ifndef (__PGCS_INTEG_METHOD__)
	#error "Please, set __PGCS_INTEG_METHOD__ to one defined in Lib_A_NINTEG_numerical_integration in macros list"
#endif*/
/*==== |End  | <-- Секция - Выбор метода численного интегрирования ===========*/

/*-------------------------------------------------------------------------*//**
 * @brief  Константа, определяющая длину большой полуоси Земли для операции
 *  обратной проекции
 */
#define PGCS_RE		(__PGCS_FPT__)6378137.0

/*-------------------------------------------------------------------------*//**
 * @brief  Переопределение числа Пи из стандартной библиотеки math.h
 */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define PGCS_PI 	(__PGCS_FPT__)M_PI

/*-------------------------------------------------------------------------*//**
 * @brief  Перечисляемый тип, задающий расположение параметров системы в векторе
 *         пространства состояний
 */
typedef enum
{
	PGCS_POS_X = 0u,		/*!< Оценка пути по оси "X" */
	PGCS_POS_Y,				/*!< Оценка пути по оси "Y" */
	PGCS_POS_Z,				/*!< Оценка пути по оси "Z" */
	/**
	 * @brief Длина вектора пространства состояний
	 */
	PGCS_LEN_STATE,
} pgcs_state_param_position_e;

#define  pgcs_fnc_status_e 					ukfmo_fnc_status_e

/*-------------------------------------------------------------------------*//**
 * @brief    Количество строк матрицы сигма-точек
 */
#define PGCS_LEN_SIGMA_ROW 					(PGCS_LEN_STATE)
/*-------------------------------------------------------------------------*//**
 * @brief    Количество столбцов матрицы сигма-точек
 */
#define PGCS_LEN_SIGMA_COL 					((PGCS_LEN_STATE * 2u) + 1u)


/*-------------------------------------------------------------------------*//**
 * @brief    Длина вектора весовых коэффициентов "среднего"
 */
#define PGCS_LEN_VECT_MEAN 					(PGCS_LEN_SIGMA_COL)
/*-------------------------------------------------------------------------*//**
 * @brief    Длина вектора весовых коэффициентов "ковариации"
 */
#define PGCS_LEN_VECT_COV					(PGCS_LEN_SIGMA_COL)


/*-------------------------------------------------------------------------*//**
 * @brief    Количество строк матриц фильтра Калмана
 */
#define PGCS_LEN_MATRIX_ROW 				PGCS_LEN_STATE
/*-------------------------------------------------------------------------*//**
 * @brief    Количество столбцов матриц фильтра Калмана
 */
#define PGCS_LEN_MATRIX_COL 				PGCS_LEN_STATE
/*#### |End  | <-- Секция - "Определение констант" ###########################*/



/*#### |Begin| --> Секция - "Определение типов" ##############################*/

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для матрицы размерностью 3x3
 */
typedef struct
{
	ukfmo_matrix_data_s mat_s;
	__PGCS_FPT__ memForMatrix[PGCS_LEN_MATRIX_ROW][PGCS_LEN_MATRIX_COL];
} pgcs_matrix_3x3_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для матрицы размерностью 3x1
 */
typedef struct
{
	ukfmo_matrix_s mat_s;
	__PGCS_FPT__ memForMatrix[PGCS_LEN_STATE][1u];
} pgcs_matrix_3_1_s;

typedef struct
{
	ukfmo_matrix_s mat_s;
	__PGCS_FPT__ memForMatrix[1u][PGCS_LEN_STATE];
} pgcs_matrix_1_3_s;

typedef struct
{
	ukfmo_matrix_s mat_s;
	__PGCS_FPT__ memForMatrix[PGCS_LEN_SIGMA_COL][1u];
} pgcs_matrix_7_1_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для матрицы размерностью 3x7
 */
typedef struct
{
	ukfmo_matrix_s mat_s;
	__PGCS_FPT__ memForMatrix[PGCS_LEN_SIGMA_ROW][PGCS_LEN_SIGMA_COL];
} pgcs_matrix_3_7_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для хранения матриц шумов "процесса" и "измерения"
 */
typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief 	Process noise covariance matrix
	 */
	pgcs_matrix_3x3_s QMat_s;

	/*------------------------------------------------------------------------*//**
	 * @brief 	Measurement noise covariance matrix
	 */
	pgcs_matrix_3x3_s RMat_s;
} pgcs_noise_matrix_s;

typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief Корень квадратный из суммы "Lаmbda" и "PGCS_LEN_STATE"
	 */
	__PGCS_FPT__ sqrtLamLen;
} pgcs_scalar_params_s;


/*-------------------------------------------------------------------------*//**
 * @brief  Структура для хранения кинематических данных
 */
typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор скорости в проекционной системе
	 */
	__PGCS_FPT__ 	flat_vel[3u];
	uint8_t			isNewVel_flag;

	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор приращения позиции в проекционной системе
	 */
	__PGCS_FPT__ 	flat_dpos[3u];

	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор начальной позиции в ДШВ-системе (Z means Zero)
	 */
	__PGCS_FPT__ 	lla_pos_zero[3u];

	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор текущей позиции в ДШВ-системе
	 */
	__PGCS_FPT__ 	lla_pos[3u];

	/*------------------------------------------------------------------------*//**
	 * @brief 	Вектор текущей позиции в ДШВ-системе, полученной с GNSS
	 */
	__PGCS_FPT__ 	lla_pos_gnss[3u];
	uint8_t			isNewPos_flag;

	/*------------------------------------------------------------------------*//**
	 * @brief 	Структура с интегралом координат в проекционной системе
	 */
	ninteg_trapz_s 	flat_pos_integ[3u];

	__PGCS_FPT__ re_c;

} pgcs_kin_data_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Структура для хранения данных, необходимых для UKF
 */
typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief Вектор пространства состояний
	 *
	 * @warning Удалить! @see "x_predict_temp_s"
	 */
	pgcs_matrix_3x3_s 	stateMat_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрицы шумов
	 */
	pgcs_noise_matrix_s noiseMatrix_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрицы ковариации ("P_k-1")
	 */
//	pgcs_matrix_3x3_s 	covMat_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Корень квадратный от матрицы ковариации ("sqrt P")
	 */
	pgcs_matrix_3x3_s 	sqrtP_apriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица распределения сигма-точек (chi_k-1)
	 */
	pgcs_matrix_3_7_s 	chiSigmaMat_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица сигма-точек (после функции преобразования) (chi_k|k-1)
	 */
	pgcs_matrix_3_7_s 	chiSigmaPostMat_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Структура, содержащая скалярные параметры фильтра
	 */
	pgcs_scalar_params_s scalar_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица x_k|k-1
	 */
	pgcs_matrix_3_1_s x_apriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица x_k
	 */
	pgcs_matrix_3_1_s x_posteriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица для хранения временных данных (используется на шаге Calculate the sigma-points)
	 */
	pgcs_matrix_3x3_s x_predict_temp_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Единичная матрица (используется на шаге Calculate the sigma-points)
	 */
	pgcs_matrix_1_3_s x_predict_temp_ones_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец весовых коэффициентов
	 */
	pgcs_matrix_7_1_s muMean_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец весовых коэффициентов
	 */
	pgcs_matrix_7_1_s muCovar_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец для хранения промежуточных результатов
	 *         (используется на шаге Calculate covariance of predicted state)
	 */
	pgcs_matrix_3_1_s chi_apriory_minus_x_apriory_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Вектор-столбец для хранения промежуточных результатов
	 *         (используется на шаге Calculate covariance of predicted state)
	 */
	pgcs_matrix_3_1_s chi_apriory_minus_x_apriory_Transpose_s;

	/*------------------------------------------------------------------------*//**
	* @brief  Матрицы для хранения промежуточных результатов
	*         (используется на шаге Calculate covariance of predicted state)
	*/
	pgcs_matrix_3x3_s resultOfMult2Matrix_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица P_k|k-1
	 */
	pgcs_matrix_3x3_s P_apriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица P_k-1
	 */
	pgcs_matrix_3x3_s P_predict_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица psi_k|k-1
	 */
	pgcs_matrix_3_7_s psi_apriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица y_k|k-1
	 */
	pgcs_matrix_3_1_s y_apriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица Pyy
	 */
	pgcs_matrix_3x3_s Pyy_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица PyyTmp
	 */
	pgcs_matrix_3x3_s PyyTmp_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица для хранения обратной от Pyy матрицы
	 */
	pgcs_matrix_3x3_s PyyInv_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица Pyy
	 */
	pgcs_matrix_3x3_s Pxy_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица коэффициентов усиления фильтра
	 */
	pgcs_matrix_3x3_s K_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Матрица для хранения обратной от матрицы коэффициентов усиления
	 */
	pgcs_matrix_3x3_s K_Transpose_s;

	/*------------------------------------------------------------------------*//**
	 * @brief  Матрица для хранения промежуточных данных (используется на шаге Calculate covariance of predicted output)
	 */
	pgcs_matrix_1_3_s psi_priory_MINUS_y_priory_TRANSPOSE;

	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец измерений, полученных от GNSS модуля
	 */
	pgcs_matrix_3_1_s y_posteriori_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Вектор-столбец "инновации"
	 */
	pgcs_matrix_3_1_s innovation_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Структура для инициализации указателей на матрицы (дли
	 *        использования библиотеки "UKFSIF")
	 */
	ukfsif_all_data_s 		ukfsifMatrixPointers_s;
} pgcs_ukf_data_s;

/*-------------------------------------------------------------------------*//**
 * @brief  Основная рабочая структура модуля PGCS, содержащая данные по
 * кинематике и UKF
 */
typedef struct
{
	pgcs_kin_data_s 		kinData_s;
	pgcs_ukf_data_s 		ukfData_s;

} pgcs_data_s;

typedef struct
{
	/*------------------------------------------------------------------------*//**
	 * @brief Скалярные параметры UKF
	 */
	ukfsif_scaling_param_s 	pgcs_scalParams_s;

	/*------------------------------------------------------------------------*//**
	 * @brief Период интегрирования измерений скорости
	 */
	__PGCS_FPT__ dt;

	/*-------------------------------------------------------------------------*//**
	 * @brief  Значения этой матрицы на шаге инициализации будут записаны в
	 *         диагональ матрицы шумов "Q"
	 */
	__PGCS_FPT__ Q_mat_a[PGCS_LEN_STATE];

	/*-------------------------------------------------------------------------*//**
	 * @brief  Значения этой матрицы на шаге инициализации будут записаны в
	 *         диагональ матрицы шумов "Q"
	 */
	__PGCS_FPT__ R_mat_a[PGCS_LEN_STATE];

	/*------------------------------------------------------------------------*//**
	 * @brief  Начальное значение вектора пространства состояний
	 */
	__PGCS_FPT__ state_a[PGCS_LEN_STATE];

} pgcs_data_init_s;


/*#### |End  | <-- Секция - "Определение типов" ##############################*/


/*#### |Begin| --> Секция - "Определение глобальных переменных" ##############*/

/*------------------------------------------------------------------------*//**
* @brief Опорный параметр для перехода из плоской системы координат в эллипсоидную
* (WGS84). Строится на основе широты стартовой точки
*/
static __PGCS_FPT__ re_c;

/*#### |End  | <-- Секция - "Определение глобальных переменных" ##############*/


/*#### |Begin| --> Секция - "Прототипы глобальных функций" ###################*/

/*------------------------------------------------------------------------*//**
* @brief Флаги наличия новых данных о скорости
*/
#define __PGCS_SetFlagVelDataUpdate() 		(pData_s->kinData_s.isNewVel_flag 	= 1u)
#define __PGCS_ReSetFlagVelDataUpdate() 	(pData_s->kinData_s.isNewVel_flag 	= 0u)
#define __PGCS_IsFlagVelDataUpdateSet() 	(pData_s->kinData_s.isNewVel_flag 	== 1u)

/*------------------------------------------------------------------------*//**
* @brief Флаги наличия новых данных о координате ДШВ, полученных с GNSS
*/
#define __PGCS_SetFlagPosDataUpdate() 		(pData_s->kinData_s.isNewPos_flag 	= 1u)
#define __PGCS_ReSetFlagPosDataUpdate() 	(pData_s->kinData_s.isNewPos_flag 	= 0u)
#define __PGCS_IsFlagPosDataUpdateSet() 	(pData_s->kinData_s.isNewPos_flag 	== 1u)


extern void
PGCS_StructInit(
    pgcs_data_init_s *pInit_s);

extern void
PGCS_Init_All(
    pgcs_data_s *pData_s,
    pgcs_data_init_s *pInit_s);

extern pgcs_fnc_status_e
PGCS_UpdatePosState(
    pgcs_data_s *pData_s);

extern void
PGCS_SetCurrentFlatVelocity(
    pgcs_data_s *pData_s,
    __PGCS_FPT__ *pVel);

extern void
PGCS_SetCurrentLLAPos(
    pgcs_data_s *pData_s,
    __PGCS_FPT__ *pLatLonAlt);

extern void
PGCS_SetZeroLLAPos(
    pgcs_data_s *pData_s,
    __PGCS_FPT__ *pLatLonAlt);

extern void
PGCS_GetProcessedLLAPos(
    pgcs_data_s *pData_s,
    __PGCS_FPT__ *pLatLonAlt);
/*#### |End  | <-- Секция - "Прототипы глобальных функций" ###################*/


/*#### |Begin| --> Секция - "Определение макросов" ###########################*/
#if defined (__UKFMO_CHEKING_ENABLE__)

/*-------------------------------------------------------------------------*//**
 * @author    Mickle Isaev
 * @date      22-авг-2019
 *
 * @brief    Макрос проверяет валидность структуры матрицы, если матрица
 *           не валидна, то макрос зацикливает программу
 *
 * @param[in]	x: 	Указатель на структуру матрицы
 *
 * @return   None
 */
__PGCS_ALWAYS_INLINE ukfmo_matrix_s*
__PGCS_CheckMatrixStructValidation(
    ukfmo_matrix_s *pData)
{
	/* Вызов макроса для проверки параметров структуры */
	__UKFMO_CheckMatrixStructValidationGeneric(
	    pData,
	    (PGCS_LEN_SIGMA_COL),
	    (PGCS_LEN_SIGMA_COL));

	return (pData);
}
#else

#define __PGCS_CheckMatrixStructValidation(x) 	(x)
#endif
/*#### |End  | <-- Секция - "Определение макросов" ###########################*/


/*#### |Begin| --> Секция - "Include - подмодули" ############################*/
/*#### |End  | <-- Секция - "Include - подмодули" ############################*/

#endif  /* LIB_A_PGCS */

/*############################################################################*/
/*################################ END OF FILE ###############################*/
/*############################################################################*/
