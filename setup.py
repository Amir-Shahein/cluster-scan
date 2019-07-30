from setuptools import setup
import pwm_scan

setup(
	name='pwm_scan',
    version=pwm_scan.__version__,
    description='Position-weight-matrix (PWM) scan through genomic sequence',
    url='https://github.com/linyc74/pwm_scan',
    author='Yu-Cheng Lin',
    author_email='linyc74@gmail.com',
    license='MIT',
    packages=['pwm_scan'],
    python_requires='>3.6',
    install_requires=['numpy', 'pandas'],
    zip_safe=False
)
