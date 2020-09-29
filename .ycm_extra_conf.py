import os

homedir = os.getenv("HOME")
def Settings( **kwargs ):
  return {
    'flags': ['-x', 'c++', 'std=c++17', '-Wall', '-Wextra', '-Wno-unused-parameter',
      '-I', '{}/.local/include'.format(homedir), '-I', 'include', '-I', '../include' , '-lIOH']
  }
