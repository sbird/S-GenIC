#S-GenIC is Simeon Bird's modifications to N-GenIC, in C++.
project('S-GenIC', 'cpp')

#External dependencies
cpp = meson.get_compiler('cpp')
if cpp.has_argument('-fopenmp')
    add_project_arguments('-fopenmp', language: 'cpp')
    omp_dep = cpp.find_library('gomp', required : false)
endif
add_project_arguments('-DHAVE_BGFL', language: 'cpp')
fftw3th = cpp.find_library('fftw3f_threads', required : false)
fftw3fth = cpp.find_library('fftw3_threads', required : false)
fftw3f = dependency('fftw3f')
fftw3 = dependency('fftw3')
gread = subproject('GadgetReader')
deps = [fftw3, fftw3f, fftw3th, fftw3fth, dependency('gsl'),gread.get_variable('gwrite_dep'),omp_dep]

src = ['cosmology.cpp', 'read_param.cpp',  'thermalvel.cpp', 'displacement.cpp', 'part_data.cpp',  'save.cpp', 'gsl_spline_wrapper.cpp',  'power.cpp']
executable('N-GenIC', sources: ['main.cpp']+src, dependencies: deps)
testdep = [dependency('boost', modules: 'test'),]
testexe = executable('btest', sources: ['test.cpp']+src, dependencies: deps + testdep)
test('Tests',testexe, workdir: meson.source_root())
