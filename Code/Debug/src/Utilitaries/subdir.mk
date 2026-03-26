################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Utilitaries/FluxPoint2D.cpp \
../src/Utilitaries/Laplace_Stehfest.cpp \
../src/Utilitaries/LinearSystem.cpp \
../src/Utilitaries/Point_Cgal.cpp \
../src/Utilitaries/RandomNumber.cpp \
../src/Utilitaries/Segment.cpp \
../src/Utilitaries/Structures.cpp \
../src/Utilitaries/UblasStructures.cpp \
../src/Utilitaries/scale.cpp 

OBJS += \
./src/Utilitaries/FluxPoint2D.o \
./src/Utilitaries/Laplace_Stehfest.o \
./src/Utilitaries/LinearSystem.o \
./src/Utilitaries/Point_Cgal.o \
./src/Utilitaries/RandomNumber.o \
./src/Utilitaries/Segment.o \
./src/Utilitaries/Structures.o \
./src/Utilitaries/UblasStructures.o \
./src/Utilitaries/scale.o 

CPP_DEPS += \
./src/Utilitaries/FluxPoint2D.d \
./src/Utilitaries/Laplace_Stehfest.d \
./src/Utilitaries/LinearSystem.d \
./src/Utilitaries/Point_Cgal.d \
./src/Utilitaries/RandomNumber.d \
./src/Utilitaries/Segment.d \
./src/Utilitaries/Structures.d \
./src/Utilitaries/UblasStructures.d \
./src/Utilitaries/scale.d 


# Each subdirectory must supply rules for building sources it contributes
src/Utilitaries/%.o: ../src/Utilitaries/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include -I/usr/include/c++/5/bits -I/home/marwa/libraries/boost_1_61_0 -I/home/marwa/libraries/boost-numeric-bindings -I/home/marwa/libraries/SuiteSparse/UMFPACK/Include -I/home/marwa/libraries/SuiteSparse/AMD/Include -I/home/marwa/libraries/CGAL-4.12-beta1/include -I/home/marwa/libraries/SuiteSparse/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


