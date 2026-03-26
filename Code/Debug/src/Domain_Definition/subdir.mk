################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Domain_Definition/DFNComputation.cpp \
../src/Domain_Definition/Domain.cpp \
../src/Domain_Definition/FractureMesh.cpp \
../src/Domain_Definition/HydraulicProperties.cpp \
../src/Domain_Definition/NetworkMeshes.cpp 

OBJS += \
./src/Domain_Definition/DFNComputation.o \
./src/Domain_Definition/Domain.o \
./src/Domain_Definition/FractureMesh.o \
./src/Domain_Definition/HydraulicProperties.o \
./src/Domain_Definition/NetworkMeshes.o 

CPP_DEPS += \
./src/Domain_Definition/DFNComputation.d \
./src/Domain_Definition/Domain.d \
./src/Domain_Definition/FractureMesh.d \
./src/Domain_Definition/HydraulicProperties.d \
./src/Domain_Definition/NetworkMeshes.d 


# Each subdirectory must supply rules for building sources it contributes
src/Domain_Definition/%.o: ../src/Domain_Definition/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include -I/usr/include/c++/5/bits -I/home/marwa/libraries/boost_1_61_0 -I/home/marwa/libraries/boost-numeric-bindings -I/home/marwa/libraries/SuiteSparse/UMFPACK/Include -I/home/marwa/libraries/SuiteSparse/AMD/Include -I/home/marwa/libraries/CGAL-4.12-beta1/include -I/home/marwa/libraries/SuiteSparse/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


