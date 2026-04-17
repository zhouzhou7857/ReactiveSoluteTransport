################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/PERFORM.cpp 

OBJS += \
./src/PERFORM.o 

CPP_DEPS += \
./src/PERFORM.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I /home/zhouw/softs/ -I /home/zhouw/softs/boost/py37/gcc95/1.76.0/include/ -I /home/zhouw/softs/cgal/gcc95/5.5.2/include/ -I ./../../../Library -I /home/zhouw/softs/suitesparse/gcc95/5.10.1/include/ -I /usr/include/GL/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

