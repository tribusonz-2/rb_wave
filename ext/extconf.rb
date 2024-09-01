# frozen_string_literal: true
require 'mkmf'

have_func('cyl_bessel_i0', 'math.h')

$INCFLAGS << ' -Iinclude'

create_makefile('wave')
