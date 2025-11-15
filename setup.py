from setuptools import setup, find_namespace_packages


setup(
    name='auto-cpo',
    version='0.3.1',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "auto-cpo = auto_cpo.__main__:main",
        ]
    },
    scripts=[],
    install_requires=[
        'requests~=2.32',
        'jinja2~=3.1',
        
    ],
    package_data={
        "auto_cpo": ["templates/*.html"],
    },
    include_package_data=True,
    description='Automated analysis of carbapenemase-producing organism (CPO) sequence data',
    url='https://github.com/BCCDC-PHL/auto-cpo',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',

    keywords=[],
    zip_safe=False
)
