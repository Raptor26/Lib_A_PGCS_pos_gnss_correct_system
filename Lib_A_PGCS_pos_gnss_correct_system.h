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
#include "Lib_A_UKFMO_ukf_matrix_operations.h"
#include "Lib_A_UKFSIF_ukf_standart_init_fnc.h"
/*==== |End  | <-- Секция - "Extern libraries" ===============================*/
/*#### |End  | <-- Секция - "Include" ########################################*/


/*#### |Begin| --> Секция - "Определение констант" ###########################*/
#if defined (__PGCS_EXTERN_MODE_ENABLE__)
	#include "macros_definitions.h"
#endif

/*==== |Begin| --> Секция определения типа числа с плавающей точкой ==========*/
#if !defined (__PGCS_FPT__)
	#error "Please, set __PGCS_FPT__ float or double in macros list"
#endif

#if !defined (__PGCS_FPT_SIZE__)
	#error "Please, set __PGCS_FPT_SIZE__ 4 (that mean float) or 8 (that mean double) in macros list"
#endif

#if     __PGCS_FPT_SIZE__ == 4

#elif   __PGCS_FPT_SIZE__ == 8

#else
	#error "Your compiler uses a non-standard floating point size"
#endif
/*==== |End  | <-- Секция определения типа числа с плавающей точкой ==========*/

/*==== |Begin| --> Секция - Макросы для встраиваемых функций =================*/
#if defined (__GNUC__)

	/* inline*/
	#ifndef __PGCS_INLINE
		#define __PGCS_INLINE          inline
	#endif

	/* static inline */
	#ifndef __PGCS_STATIC_INLINE
		#define __PGCS_STATIC_INLINE   static inline
	#endif

	/* always inline */
	#ifndef __PGCS_ALWAYS_INLINE
		#define __PGCS_ALWAYS_INLINE    inline __attribute__((always_inline)) static
	#endif

	/* force inline */
	#ifndef __PGCS_FORCE_INLINE
		#define __PGCS_FORCE_INLINE    inline __attribute__((always_inline))
	#endif

#else
	#define __PGCS_INLINE
	#define __PGCS_STATIC_INLINE   static
	#define __PGCS_ALWAYS_INLINE
#endif
/*==== |End  | <-- Секция - Макросы для встраиваемых функций =================*/

/*==== |Begin| --> Секция - Расположение функций библиотеки в специальной
 *                          области памяти ===================================*/
#if defined (__PGCS_FNC_ONCE_MEMORY_LOCATION_NAME__)
	#if defined (__GNUC__)
		#define __PGCS_FNC_ONCE_MEMORY_LOCATION  __attribute__ ((section(__PGCS_FNC_ONCE_MEMORY_LOCATION_NAME__)))
	#else
		#error "You defined the name of the memory area for the function location, but the type of your compiler is not supported by the library. You can delete the macro definition __PGCS_FNC_ONCE_MEMORY_LOCATION_NAME__ or extend the macro definition __PGCS_FNC_ONCE_MEMORY_LOCATION for your compiler type"
	#endif
#else
	#define __PGCS_FNC_ONCE_MEMORY_LOCATION
#endif

#if defined (__PGCS_FNC_LOOP_MEMORY_LOCATION_NAME__)
	#if defined (__GNUC__)
		#define __PGCS_FNC_LOOP_MEMORY_LOCATION  __attribute__ ((section(__PGCS_FNC_LOOP_MEMORY_LOCATION_NAME__)))
	#else
		#error "You defined the name of the memory area for the function location, but the type of your compiler is not supported by the library. You can delete the macro definition __PGCS_FNC_LOOP_MEMORY_LOCATION_NAME__ or extend the macro definition __PGCS_FNC_LOOP_MEMORY_LOCATION for your compiler type"
	#endif
#else
	#define __PGCS_FNC_LOOP_MEMORY_LOCATION
#endif
/*==== |End  | <-- Секция - Расположение функций библиотеки в специальной
 *                          области памяти ===================================*/
/*#### |End  | <-- Секция - "Определение констант" ###########################*/


/*#### |Begin| --> Секция - "Определение типов" ##############################*/
typedef struct
{

} pgcs_data_s;

typedef struct
{

} pgcs_data_init_s;
/*#### |End  | <-- Секция - "Определение типов" ##############################*/


/*#### |Begin| --> Секция - "Определение глобальных переменных" ##############*/
/*#### |End  | <-- Секция - "Определение глобальных переменных" ##############*/


/*#### |Begin| --> Секция - "Прототипы глобальных функций" ###################*/
extern void
PGCS_StructInit(
  pgcs_data_init_s *pInit_s);

extern void
PGCS_Init_All(
  pgcs_data_s *pData_s,
  pgcs_data_init_s *pInit_s);

extern void
PGCS_UpdatePosState(
  pgcs_data_s *pData_s);
extern void
PGCS_CopyVelInWorldFrame(
  pgcs_data_s *pData_s,
  __PGCS_FPT__ *pVel);

extern void
PGCS_CopyLatLOnAltInWorldFrame(
  pgcs_data_s *pData_s,
  __PGCS_FPT__ *pLatLonAlt);
/*#### |End  | <-- Секция - "Прототипы глобальных функций" ###################*/


/*#### |Begin| --> Секция - "Определение макросов" ###########################*/
/*#### |End  | <-- Секция - "Определение макросов" ###########################*/


/*#### |Begin| --> Секция - "Include - подмодули" ############################*/
/*#### |End  | <-- Секция - "Include - подмодули" ############################*/

#endif  /* LIB_A_PGCS */

/*############################################################################*/
/*################################ END OF FILE ###############################*/
/*############################################################################*/
