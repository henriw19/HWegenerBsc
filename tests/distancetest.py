from ldpc.code_util import compute_code_distance

from hwbsc.color_code import hexa_color_pcm_generator

pcm = hexa_color_pcm_generator(12,6)[0]
d = compute_code_distance(pcm)
print(d)