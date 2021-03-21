# rebuild pyf file
# f2py xztrin21.f -m opacities -h opacities.pyf
# build C extension
#python3 $(which f2py) opacities.pyf
# after that, copy opacitiesmodule.c to ..
# build directly
python3 $(which f2py) -c opacities.pyf xztrin21.f
