################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Chemistry/Chemistry.cpp \

OBJS += \
./src/Chemistry/Chemistry.o \

CPP_DEPS += \
./src/Chemistry/Chemistry.d \


# Each subdirectory must supply rules for building sources it contributes
src/Chemistry/%.o: ../src/Chemistry/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I /home/zhouw/softs/ -I /home/zhouw/softs/boost/py37/gcc95/1.76.0/include/ -I /home/zhouw/softs/cgal/gcc95/5.5.2/include/ -I ./../../../Library -I /usr/include/suitesparse/ -I /usr/include/GL/ -frounding-math -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


