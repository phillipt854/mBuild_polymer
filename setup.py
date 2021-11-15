from setuptools import setup

setup(
      name="mbuild_polymer",
      install_requires="mbuild",
      entry_points={
                    "mbuild.plugins":[ "polymer = mbuild_polymer.mbuild_polymer:polymer_box",
                                       "general_polymer = mbuild_polymer.mbuild_general_polymer:general_polymer"
                        ]
                    },
                    py_modules=["mbuild_polymer"],
                        )
