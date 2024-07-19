# frozen_string_literal: true
require 'mkmf'

have_func('cyl_bessel_i0', 'math.h')

create_makefile('wave')
