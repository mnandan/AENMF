################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../deriveAE.o \
../fileInt.o \
../getFactors.o 

CPP_SRCS += \
../AENMF.cpp \
../deriveAE.cpp \
../fileInt.cpp \
../getFactors.cpp 

OBJS += \
./AENMF.o \
./deriveAE.o \
./fileInt.o \
./getFactors.o 

CPP_DEPS += \
./AENMF.d \
./deriveAE.d \
./fileInt.d \
./getFactors.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


