#!/usr/bin/env python

from __future__ import print_function

import os, sys

if len(sys.argv) < 2:
    print('USAGE: {} path_to_normalize'.format(sys.argv[0]))
    sys.exit(1)

script_name            = os.path.basename(sys.argv[0])
normalization_function = None

if script_name.startswith('abs_'):
    normalization_function = os.path.abspath
elif script_name.startswith('real_path'):
    normalization_function = os.path.realpath
else:
    print('Unknown function', file = sys.stderr)
    sys.exit(1)

print(normalization_function(sys.argv[1]))
