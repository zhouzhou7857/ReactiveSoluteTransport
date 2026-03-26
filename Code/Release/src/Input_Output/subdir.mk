################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Input_Output/BoundaryConditions.cpp \
../src/Input_Output/Parameters.cpp \
../src/Input_Output/Results.cpp 

OBJS += \
./src/Input_Output/BoundaryConditions.o \
./src/Input_Output/Parameters.o \
./src/Input_Output/Results.o 

CPP_DEPS += \
./src/Input_Output/BoundaryConditions.d \
./src/Input_Output/Parameters.d \
./src/Input_Output/Results.d 


# Each subdirectory must supply rules for building sources it contributes
src/Input_Output/%.o: ../src/Input_Output/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I /home/zhouw/softs/ -I /home/zhouw/softs/boost/py37/gcc95/1.76.0/include/ -I /home/zhouw/softs/cgal/gcc95/5.5.2/include/ -I ./../../../Library -I /usr/include/suitesparse/ -I /usr/include/GL/ -frounding-math -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


