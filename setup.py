from setuptools import setup

setup(
      name="mbuild_polymer",
      install_requires="mbuild",
      entry_points={
                    "mbuild.plugins":[ "polymer = mbuild_polymer.mbuild_polymer:polymer_box"
                        ]
                    },
                    py_modules=["mbuild_polymer"],
                        )
