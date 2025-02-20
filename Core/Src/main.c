/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2024 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "adc.h"
#include "dma.h"
#include "i2c.h"
#include "tim.h"
#include "usart.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "stdio.h"
#include "OLED.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
uint32_t timeMark = 0;//全局时间量，spend 10us add 1
uint16_t Voltage = 0;
uint16_t Current = 0;
uint16_t Power = 0;
uint16_t Freq = 0;
uint16_t voltageAndCurrentValue[2048] = {0};

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
//清空数组
void clearArray(uint16_t* Array,uint32_t Length,uint16_t defaultValue){
  for (int i = 0; i < Length; i++)
  {
    *Array = defaultValue;
    Array++;
  }
  
}
//串口重定�?
int fputc(int ch, FILE *f){
  HAL_UART_Transmit(&huart1, (uint8_t *)&ch, 1, 0xffff);
  return ch;
}
//输入1个十进制的数，返回十进制数的长度
uint8_t numBitCount(uint32_t number){
  uint8_t Length = 0;
  if (number == 0) return 1;  
  while (number)
  {
    number = number/10;
    Length++;
  }
  return Length;
}
//更新数据
void updateData(void){
  OLED_Clear();
  OLED_ShowString(0,0,"Voltage:",OLED_8X16);
  OLED_ShowString(0,16,"Current:",OLED_8X16);
  OLED_ShowString(0,32,"Power:",OLED_8X16);
  OLED_ShowString(0,48,"Freq:",OLED_8X16);
  OLED_ShowNum(64,0,Voltage,numBitCount(Voltage),OLED_8X16);
  OLED_ShowNum(64,16,Current,numBitCount(Current),OLED_8X16);
  OLED_ShowNum(64,32,Power,numBitCount(Power),OLED_8X16);
  OLED_ShowNum(64,48,Freq,numBitCount(Freq),OLED_8X16);
  OLED_ShowString(112,0,"V",OLED_8X16);
  OLED_ShowString(112,16,"mA",OLED_8X16);
  OLED_ShowString(112,32,"W",OLED_8X16);
  OLED_ShowString(112,48,"Hz",OLED_8X16);
  OLED_Update();
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_ADC1_Init();
  MX_I2C1_Init();
  MX_USART1_UART_Init();
  MX_TIM14_Init();
  /* USER CODE BEGIN 2 */
  OLED_Init();
  HAL_TIM_Base_Start_IT(&htim14);
  HAL_ADCEx_Calibration_Start(&hadc1);
  HAL_ADC_Start_DMA(&hadc1,(uint32_t*)voltageAndCurrentValue,sizeof(voltageAndCurrentValue)/sizeof(uint16_t));
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  uint16_t voltageMax = 0;
  uint16_t currentMax = 0;
  uint16_t maxVoltageLocation = 0;
  uint16_t minVoltageLocation = 0;

  double freqAccumulator = 0;
  double currentAccumulator = 0;
  double voltageAccumulator = 0;
  uint16_t AccumulateCount = 1;

  while (1)
  {
    HAL_GPIO_WritePin(LED1_GPIO_Port,LED1_Pin,GPIO_PIN_RESET);//正常工作,常亮

    if (voltageAndCurrentValue[sizeof(voltageAndCurrentValue)/sizeof(uint16_t) - 1] != 0x0000)
    {
      HAL_ADC_Stop_DMA(&hadc1);
      //电压
      for (int i = 0; i < sizeof(voltageAndCurrentValue)/sizeof(voltageAndCurrentValue[0]); i += 2)
      {
        if (voltageAndCurrentValue[i] > voltageMax)
        {
          voltageMax = voltageAndCurrentValue[i];
        }
      }
      voltageAccumulator += ((voltageMax* 600)/4096)-240;
      Voltage = voltageAccumulator/(1.414 * AccumulateCount);
      //电流
      for (int i = 1; i < sizeof(voltageAndCurrentValue)/sizeof(voltageAndCurrentValue[0]); i += 2)
      {
       if (voltageAndCurrentValue[i] > currentMax)
       {
         currentMax = voltageAndCurrentValue[i];
       }
      }
      Current = (currentMax * 3300/4096) - 2243;
			if(Current > 3000) Current = 0;
      currentAccumulator += Current;
			Current = currentAccumulator/(1.414 * AccumulateCount);
      Power = Voltage * Current / 1000; 
      //频率
			uint16_t max = 0;
			uint16_t min = voltageAndCurrentValue[0];
      for (int i = 0; i < sizeof(voltageAndCurrentValue)/sizeof(voltageAndCurrentValue[0]); i += 2)
      {
        if (voltageAndCurrentValue[i] > max)
        {
          max =  voltageAndCurrentValue[i];
          maxVoltageLocation = i;
        }
        if (voltageAndCurrentValue[i] < min)
        {
          min = voltageAndCurrentValue[i];
          minVoltageLocation = i;
        }
      }
      uint16_t distance = (maxVoltageLocation > minVoltageLocation)?(maxVoltageLocation - minVoltageLocation):(minVoltageLocation - maxVoltageLocation);
      double Frequeny = 1000000/((distance - 1) * 10.8125);
      freqAccumulator += Frequeny;
      Freq = freqAccumulator / (AccumulateCount * 2);
			AccumulateCount++;

      //累加值清零,防止溢出
      if(AccumulateCount == 2000){
        AccumulateCount = 1;
        freqAccumulator = 0;
        currentAccumulator = 0;
        voltageAccumulator = 0;
      }
      // for (int i = 0; i < sizeof(voltageAndCurrentValue)/sizeof(voltageAndCurrentValue[0]); i += 2){
      //   printf("%d ",voltageAndCurrentValue[i]);
      // }
      // HAL_Delay(5000);
      clearArray(voltageAndCurrentValue,sizeof(voltageAndCurrentValue)/sizeof(voltageAndCurrentValue[0]),0x0000);
      HAL_ADC_Start_DMA(&hadc1,(uint32_t*)voltageAndCurrentValue,sizeof(voltageAndCurrentValue)/sizeof(uint16_t));
    }
    updateData();//刷新数据到OLED
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  HAL_PWREx_ControlVoltageScaling(PWR_REGULATOR_VOLTAGE_SCALE1);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = RCC_PLLM_DIV1;
  RCC_OscInitStruct.PLL.PLLN = 16;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLR = RCC_PLLR_DIV2;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
  if (htim == &htim14)
  {
    timeMark++;
  }
}
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
    //出现硬件上错误，闪烁
    HAL_GPIO_WritePin(LED1_GPIO_Port,LED1_Pin,GPIO_PIN_RESET);
    HAL_Delay(200);
    HAL_GPIO_WritePin(LED1_GPIO_Port,LED1_Pin,GPIO_PIN_SET);
    HAL_Delay(200);
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
