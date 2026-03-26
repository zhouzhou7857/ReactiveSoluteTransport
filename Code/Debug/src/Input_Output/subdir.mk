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
	g++ -I/usr/include -I/usr/include/c++/5/bits -I/home/marwa/libraries/boost_1_61_0 -I/home/marwa/libraries/boost-numeric-bindings -I/home/marwa/libraries/SuiteSparse/UMFPACK/Include -I/home/marwa/libraries/SuiteSparse/AMD/Include -I/home/marwa/libraries/CGAL-4.12-beta1/include -I/home/marwa/libraries/SuiteSparse/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


