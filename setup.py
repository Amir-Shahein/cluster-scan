from setuptools import setup

setup(name='pwm_scan',
      version='1.3',
      description='Position-weight-matrix (PWM) scan through genomic sequence',
      url='https://github.com/linyc74/pwm_scan',
      author='Yu-Cheng Lin',
      author_email='linyc74@gmail.com',
      license='MIT',
      packages=['pwm_scan'],
      install_requires=['numpy', 'pandas'],
      zip_safe=False)
