import itertools
import math
from matplotlib import rcParams
import scipy.optimize
import numpy as np
import scipy.integrate
import cubature
import scipy.optimize
from scipy.optimize import minimize

rcParams.update({'figure.autolayout': True})


def integrate(func, start, stop, use_romberg=False):
	"""
	Integrates the given function between the start and stop points.

	Parameters:
	- func (callable): The function to integrate.
	- start (float): The starting point of integration.
	- stop (float): The stopping point of integration.
	- use_romberg (bool, optional): If True, use romberg method. If False, use scipy's quad method. Default is False.

	Returns:
	- float: The result of the integration.
	"""
	if use_romberg:
		# Using romberg integration method with specified error tolerances and maximum evaluations
		#res = cubature.cubature(func, 1, 1, [start], [stop], abserr=1e-9, relerr=1e-9, maxEval=3e4)[0][0]
		#res = scipy.integrate.romberg(func, start, stop, tol=1e-15, rtol=1e-15, divmax=5)
		res = scipy.integrate.quad(func, start, stop, epsabs=1e-9, epsrel=1e-9, limit=100)[0]
	else:
		# Using scipy's quad integration method with specified error tolerances and limit
		res = scipy.integrate.quad(func, start, stop, epsabs=1e-12, epsrel=1e-12, limit=100)[0]
	return res


class Repair:
	"""
	Represents a Repair variable set to a certain value.
	"""
	
	def __init__(self, val):
		"""
		Initialize the Repair object with a given value, representing the current repair investment.

		Parameters:
		- val (float): The current repair investment (between 0 and rmax).
		"""
		self.val = val
	
	def ret_val(self):
		"""
		Returns the current repair investment.

		Returns:
		- float: The repair investment.
		"""
		return self.val


class InverseTOff:
	"""
	A trade-off function for the rate of aging that takes the form beta*((rmax/r)-1)
	"""
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta

	def ret_val(self):
		"""
		Calculates the rate of aging based on the current repair investment.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * ((self.rmax / self.repair.ret_val()) - 1)


class AsymTOff:
	def __init__(self, repair, rmax, beta):
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * ((self.rmax - self.repair.ret_val()) ** 4)


class GompertzRepAMR_Asym:
	def __init__(self, A, toff):
		self.A = A
		self.toff = toff
	
	def get_IMR(self):
		return self.A
	
	def get_AMR(self):
		#return self.toff.ret_val()
		return 0.2*(1 - math.exp(-self.toff.ret_val()))

	def get_scaled_AMR(self):
		return 0.2*(1 - math.exp(-self.toff.ret_val()))

	def ret_val_at_time(self, x):
		return self.get_IMR() * math.exp(self.get_scaled_AMR() * x)
	
	def integral_form(self, x, a, decay_start=0):
		if self.get_scaled_AMR() == 0 or self.get_IMR() == 0:
			return math.exp(-self.get_IMR() * (x - a))
		else:
			if self.get_scaled_AMR() * (x - decay_start) > 100:
				return 0
			else:
				return math.exp(-self.get_IMR() * (math.exp(self.get_scaled_AMR() * (x - decay_start)) -
				                                   math.exp(self.get_scaled_AMR() * (a - decay_start))) / self.get_AMR())


class InverseTOffPlei:
	"""
	A trade-off function for the rate of aging that takes the form beta*((rmax/(r+plei))-1)
	
	Same as the InverseTOff function with an additional plei value added to r, representing some level of free
	pleiotropic repair that is added on top of the invested repair.
	"""
	def __init__(self, repair, rmax, beta, plei):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		- plei (float): Pleiotropic repair level.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
		self.plei = plei

	def ret_val(self):
		"""
		Compute the inverse trade-off value with respect to the repair value,
		the additional 'plei' value, and the maximum bounding value.

		Returns:
		float: The computed inverse trade-off value.
		"""
		if self.repair.ret_val() + self.plei >= self.rmax:
			return 0
		return self.beta * ((self.rmax / (self.repair.ret_val() + self.plei)) - 1)


class InverseTOffMinus:
	"""
	A trade-off function for the rate of aging that takes the form beta*((1/r)-(1/rmax))
	"""
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta

	def ret_val(self):
		"""
		Compute the inverse trade-off difference between the repair value
		and the maximum bounding value.

		Returns:
		float: The computed inverse trade-off difference.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return (1 / self.repair.ret_val()) - (1 / self.rmax)

class SqrtTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (1 - sqrt(repair + 1 - rmax))
	"""
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta

	def ret_val(self):
		"""
		Calculate the rate of aging based on the current repair investment using
		a square root function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (1 - (self.repair.ret_val() + 1 - self.rmax) ** 0.5)


class SquareTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (1 - (repair + 1 - rmax)^2)
	"""
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta

	def ret_val(self):
		"""
		Calculate the rate of aging based on the current repair investment using
		a squared function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (1 - (self.repair.ret_val() + 1 - self.rmax) ** 2)


class SharpTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (1 - (repair^20 / rmax^20))
	"""
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta

	def ret_val(self):
		"""
		Calculate the rate of aging based on the current repair investment using
		a sharply rising power function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (1 - (self.repair.ret_val() ** 20 / self.rmax ** 20))


class ScalingTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (rmax^5 - repair^5).
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (self.rmax ** 5 - self.repair.ret_val() ** 5)


class ExpTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (exp(-(repair - rmax)) - 1).
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (math.exp(-(self.repair.ret_val() - self.rmax)) - 1)


class LinearTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (rmax - repair).
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (self.rmax - self.repair.ret_val())


class LinearMinusTOff:
	"""
	A trade-off function for the rate of aging that takes the form
	beta * (1 - (repair / rmax)).
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value for which the rate of aging is zero.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		return self.beta * (1 - (self.repair.ret_val() / self.rmax))


class LinearDiscontinuousTOff:
	"""
	A discontinuous trade-off function for the rate of aging.
	- Returns beta * (rmax - 2 * repair) for repair < rmax/2
	- Returns beta * (rmax/2 - repair/2) for rmax/2 <= repair < rmax
	- Returns 0 for repair >= rmax
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value which affects the computation of the trade-off.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		elif self.repair.ret_val() >= (self.rmax / 2):
			return self.beta * ((self.rmax / 2) - (self.repair.ret_val() / 2))
		else:
			return self.beta * (self.rmax - (2 * self.repair.ret_val()))


class BinaryTOff:
	"""
	A trade-off function for the rate of aging that produces a binary output:
	- Returns beta if repair < rmax
	- Returns 0 for repair >= rmax
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair threshold value which determines the trade-off output.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		else:
			return self.beta


class StepTOff:
	"""
	A trade-off function for the rate of aging that produces step-wise decrements in output
	as repair approaches rmax. Takes the same form as InverseTOff
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value that serves as the denominator in the step function.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		else:
			return self.beta * ((1 / ((math.floor(10 * self.repair.ret_val() / self.rmax) / 10) + 0.01)) - 1)


class SmoothStepTOff:
	"""
	A trade-off function for the rate of aging that produces step-wise decrements in output
	as repair approaches rmax. Takes the same form as InverseTOff
	"""
	
	def __init__(self, repair, rmax, beta):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- rmax (float): Repair value that serves as the denominator in the step function.
		- beta (float): Scaling factor for the rate of aging.
		"""
		self.repair = repair
		self.rmax = rmax
		self.beta = beta
	
	def ret_val(self):
		"""
		Calculate the rate of aging using the specified trade-off function.

		Returns:
		float: The rate of aging.
		"""
		if self.repair.ret_val() >= self.rmax:
			return 0
		else:
			factor = self.rmax / (0.1 ** (1 / 3) + 0.5)
			return self.beta * (-(((((self.repair.ret_val() / factor) - 0.5) ** 3)) / 0.1)+1)


class ReproductiveScaling:
	"""
	The reproductive side of the trade-off function that scales reproduction down by (1-r)
	"""
	def __init__(self, repair):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		"""
		self.repair = repair
	
	def ret_val(self):
		"""
		Calculate the scaling factor of reproduction based on the current repair investment

		Returns:
		float: The scaling factor
		"""
		return 1 - self.repair.ret_val()


class ReproductiveScalingExp:
	"""
	The same as ReproductiveScaling, but with reproduction scaled by e^(-r)
	
	If we have h(r)=e^(-r) as our scaling function, we could redefine our functions to be
	H(r)=(1-r), B(r) = b(h^-1(H(r)), and A(r)=a(h^-1(H(r)) without changing anything. Thus,
	
	This function (and the others) were used as sanity checks for this fact. However, our paper
	defaults to leaving h(r)=1-r and focusing on differences in b(r) and a(r).
	"""
	def __init__(self, repair):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		"""
		self.repair = repair
	
	def ret_val(self):
		"""
		Calculate the scaling factor of reproduction based on the current repair investment

		Returns:
		float: The scaling factor
		"""
		return math.exp(-self.repair.ret_val())


class ReproductiveScalingDrain:
	"""
	A ReproductiveScaling function with an energy drain added, scaling reproduction by (1-r-d)
	"""
	
	def __init__(self, repair, d):
		"""
		Init function

		Parameters:
		- repair (object): The Repair variable.
		- d (float): The additive drain to energy stores
		"""
		self.repair = repair
		self.d = d
	
	def ret_val(self):
		"""
		Calculate the scaling factor of reproduction based on the current repair investment

		Returns:
		float: The scaling factor
		"""
		return 1 - self.repair.ret_val() - self.d


class ReproductiveScalingMulti:
	"""
	A ReproductiveScaling object that takes multiple repair variables, scaling reproduction by (1-r1-r2-r3 ...)
	"""
	def __init__(self, repairs):
		"""
		Init function

		Parameters:
		- repairs (list): A list of Repair variables.
		"""
		self.repairs = repairs
	
	def ret_val(self):
		"""
		Calculate the scaling factor of reproduction based on the current repair investment

		Returns:
		float: The scaling factor
		"""
		rep_sums = sum([x.ret_val() for x in self.repairs])
		return 1 - rep_sums


class GompertzCons:
	"""
	An object that represents a Gompertz function with a constant (not influenced by repair investment)
	value for the IMR and the AMR (given by A and B respectively). The final form is A*e^(Bx).
	"""
	def __init__(self, A, B):
		"""
		Init function

		Parameters:
		- A (float): The IMR
		- B (float): The AMR
		"""
		self.A = A
		self.B = B
	
	def get_IMR(self):
		"""
		Calculate the IMR (constant here)

		Returns:
		float: The IMR
		"""
		return self.A
	
	def get_AMR(self):
		"""
		Calculate the AMR (constant here)

		Returns:
		float: The AMR
		"""
		return self.B
	
	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age
		
		Returns:
		float: The current mortality rate
		"""
		return self.get_IMR() * math.exp(self.get_AMR() * x)
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.
		
		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins
		
		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		if self.get_AMR() == 0 or self.get_IMR() == 0:
			return math.exp(-self.get_IMR() * (x - a))
		else:
			if self.get_AMR() * (x - decay_start) > 100:
				return 0
			else:
				return math.exp(-self.get_IMR() * (math.exp(self.get_AMR() * (x - decay_start)) -
													 math.exp(self.get_AMR() * (a - decay_start))) / self.get_AMR())


class GompertzRepAMR:
	"""
	An object that represents a Gompertz function with a constant (not influenced by repair investment)
	value for the IMR. The AMR is determined by the trade-off function based on repair investment.
	The final form is A*e^(b(r)*x). Note that the beta term is factored into b(r) here.
	"""
	def __init__(self, A, toff):
		"""
		Init function

		Parameters:
		- A (float): The IMR
		- toff (object): The Tradeoff function object (like InverseTOff)
		"""
		self.A = A
		self.toff = toff
	
	def get_IMR(self):
		"""
		Calculate the IMR (constant here)

		Returns:
		float: The IMR
		"""
		return self.A
	
	def get_AMR(self):
		"""
		Calculate the AMR based on the toff function

		Returns:
		float: The AMR
		"""
		return self.toff.ret_val()
	
	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age

		Returns:
		float: The current mortality rate
		"""
		return self.get_IMR() * math.exp(self.get_AMR() * x)
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		if self.get_AMR() == 0 or self.get_IMR() == 0:
			return math.exp(-self.get_IMR() * (x - a))
		else:
			if self.get_AMR() * (x-decay_start) > 100:
				return 0
			else:
				return math.exp(-self.get_IMR() * (math.exp(self.get_AMR() * (x - decay_start)) -
													 math.exp(self.get_AMR() * (a - decay_start))) / self.get_AMR())


class GompertzRepIMR:
	"""
	An object that represents a Gompertz function with a constant (not influenced by repair investment)
	value for the AMR. The IMR is determined by the trade-off function based on repair investment.
	The final form is b(r)*e^(B*x). Note that the beta term is factored into b(r) here.
	"""
	def __init__(self, toff, B):
		"""
		Init function

		Parameters:
		- toff (object): The Tradeoff function object (like InverseTOff)
		- B (float): The AMR
		"""
		self.toff = toff
		self.B = B
	
	def get_IMR(self):
		"""
		Calculate the IMR based on the toff function

		Returns:
		float: The IMR
		"""
		return self.toff.ret_val()
	
	def get_AMR(self):
		"""
		Calculate the AMR (constant here)

		Returns:
		float: The AMR
		"""
		return self.B
	
	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age

		Returns:
		float: The current mortality rate
		"""
		return self.get_IMR() * math.exp(self.get_AMR() * x)
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		if self.get_AMR() == 0 or self.get_IMR() == 0:
			return math.exp(-self.get_IMR() * (x - a))
		else:
			if self.get_AMR() * x > 500:
				return 0
			else:
				return math.exp(-self.get_IMR() * (math.exp(self.get_AMR() * (x - decay_start)) -
													 math.exp(self.get_AMR() * (a - decay_start))) / self.get_AMR())




class GompertzRepBoth:
	"""
	An object that represents a Gompertz function with both the AMR and IMR determined by trade-off
	functions based on some repair investment. The final form is b1(r)*e^(b2(r)*x).
	Note that the beta term is factored into b1(r) and b2(r) here.
	
	This class is used to mimic cases where repair decreases the IMR and the AMR of the same Gompertz function.
	It can be used to explore cases where the same repair investment decreases IMR and AMR simultaneously or
	cases where two different repair investments decrease IMR and AMR in the same individual.
	"""
	
	def __init__(self, toff1, toff2):
		"""
		Init function

		Parameters:
		- toff1 (object): The Tradeoff function object (like InverseTOff) that controls the IMR
		- toff2 (object): The Tradeoff function object (like InverseTOff) that controls the AMR
		"""
		self.toff1 = toff1
		self.toff2 = toff2
	
	def get_IMR(self):
		"""
		Calculate the IMR based on the toff function

		Returns:
		float: The IMR
		"""
		return self.toff1.ret_val()
	
	def get_AMR(self):
		"""
		Calculate the AMR based on the toff function

		Returns:
		float: The AMR
		"""
		return self.toff2.ret_val()
	
	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age

		Returns:
		float: The current mortality rate
		"""
		return self.get_AMR() * math.exp(self.get_IMR() * x)
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		if self.get_AMR() == 0 or self.get_IMR() == 0:
			return math.exp(-self.get_IMR() * (x - a))
		else:
			if self.get_AMR() * x > 500:
				return 0
			else:
				return math.exp(-self.get_IMR() * (math.exp(self.get_AMR() * (x - decay_start)) -
													 math.exp(self.get_AMR() * (a - decay_start))) / self.get_AMR())


class WeibullRepPower:
	"""
	An object that represents a Weibull function with a constant (not influenced by repair investment)
	value for the coefficient term. The power term is determined by the trade-off function based on repair investment.
	The final form is A*x^b(r). Note that the beta term is factored into b(r) here.
	"""
	def __init__(self, A, toff):
		"""
		Init function

		Parameters:
		- toff (object): The Tradeoff function object (like InverseTOff)
		- A (float): The coefficient (constant)
		"""
		self.A = A
		self.toff = toff
	
	def get_coeff(self):
		"""
		Calculate the coeff (constant here)

		Returns:
		float: The coefficient
		"""
		return self.A
	
	def get_power(self):
		"""
		Calculate the power based on the toff function

		Returns:
		float: The power
		"""
		return self.toff.ret_val()

	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age

		Returns:
		float: The current mortality rate
		"""
		return self.get_coeff() * (x**self.get_power())
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		if self.get_power() == 0 or self.get_coeff() == 0:
			return math.exp(-self.get_coeff() * (x - a))
		else:
			return math.exp(-self.get_coeff() * (((x-decay_start)**(1+self.get_power()) ) - ((a-decay_start)**(1+self.get_power()))) / (1+self.get_power()))


class WeibullRepCoeff:
	"""
	An object that represents a Weibull function with a constant (not influenced by repair investment)
	value for the power term. The coefficient is determined by the trade-off function based on repair investment.
	The final form is b(r)*x^B. Note that the beta term is factored into b(r) here.
	"""
	def __init__(self, toff, B):
		"""
		Init function

		Parameters:
		- toff (object): The Tradeoff function object (like InverseTOff)
		- B (float): The power value
		"""
		self.toff = toff
		self.B = B
	
	def get_coeff(self):
		"""
		Calculate the coeff based on the toff function

		Returns:
		float: The coefficient
		"""
		return self.toff.ret_val()
	
	def get_power(self):
		"""
		Calculate the power (constant here)

		Returns:
		float: The power
		"""
		return self.B
	
	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age

		Returns:
		float: The current mortality rate
		"""
		return self.get_coeff() * (x ** self.get_power())
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		if self.get_coeff() == 0 or self.get_power() == 0:
			return math.exp(-self.get_coeff() * (x - a))
		else:
			return math.exp(-self.get_coeff() * (
					((x - decay_start) ** (1 + self.get_power())) - ((a - decay_start) ** (1 + self.get_power()))) / (
										1 + self.get_power()))


class ArmitageDoll:
	"""
	An object that represents the Armitage-Doll tumorigenesis model. We assume that tumorigenesis always results
	in death so that the model can be used for mortality. The probability of a tumorigenic mutation is determined
	by the trade-off function based on repair investment. The size (C) and the number of steps needed for a mutation
	(r) are constant. The final form for mortality is C/(r-1)! * p^r * x^(r-1)
	
	https://en.wikipedia.org/wiki/Armitage%E2%80%93Doll_multistage_model_of_carcinogenesis
	"""
	def __init__(self, toff, C_size, r_steps):
		"""
		Init function

		Parameters:
		- toff (object): The Tradeoff function object (like InverseTOff)
		- C_size (float): The number of cells in the individual capable of tumorigenesis
		- r_steps (float): The number of mutational steps needed for tumorigenesis
		"""
		self.toff = toff
		self.C_size = C_size
		self.r_steps = r_steps
	
	def get_p(self):
		"""
		Calculate the instantaneous probability of a mutation in an oncogene, controlled by repair

		Returns:
		float: The instantaneous probability of an oncogene
		"""
		return self.toff.ret_val()

	def risk_scale(self):
		"""
		Calculate the coefficient for the mortality function (uses the number of cells, probability of an oncogenic
		mutation, and the number of oncogenic mutations needed for tumorigenesis).

		Returns:
		float: The instantaneous probability of an oncogene
		"""
		return self.C_size * (self.get_p())**(self.r_steps)  / math.factorial(self.r_steps-1)
	
	def ret_val_at_time(self, x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age

		Returns:
		float: The current mortality rate
		"""
		return self.risk_scale() * x**(self.r_steps-1)
	
	#def integral_form(self, x, a, decay_start=0):
	#	return math.exp(-integrate(lambda x : self.ret_val_at_time(x-decay_start),a,x))

	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		The 'decay_start' is the age at which mortality should start increasing. For
		all figures in the paper, we assume that aging starts at age zero. However,
		Drenos and Kirkwood assume that this starts at the age of sexual maturity.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		return math.exp(-self.risk_scale() * ((x-decay_start)**(self.r_steps) - (a-decay_start)**(self.r_steps))/(self.r_steps))


class Extrinsic:
	"""
	An object that represents an external, constant source of mortality.
	"""
	def __init__(self, gamma):
		"""
		Init function

		Parameters:
		- gamma (float): The rate of mortality
		"""
		self.gamma = gamma
	
	def ret_val_at_time(self,x):
		"""
		Calculates the mortality rate

		Parameters:
		- x (float): The current age (not used but kept for consistency)

		Returns:
		float: The current mortality rate
		"""
		return self.gamma
	
	def integral_form(self, x, a, decay_start=0):
		"""
		Calculates e to the negative integral of the mortality rate from a to x. Serves a purpose in the
		calculation of survival curves.

		Parameters:
		- x (float): The end of the integral
		- a (float): The beginning of the integral
		- decay_start (float): The age at which aging begins (not used)

		Returns:
		float: exp(-integral of u(t)dt from a to x)
		"""
		return math.exp(-self.gamma * (x - a))


class PeakReproductiveValue:
	"""
	An object that represents the highest reproductive value possible given the current investment in repair.
	"""
	def __init__(self, toff, scale):
		"""
		Init function

		Parameters:
		- toff (object): The Tradeoff function object (like ReproductiveScaling)
		- scale (float): The scale of reproduction, such that the absolute highest value of reproduction theoretically
		                 possible with no repair is equal to the scale variable.
		"""
		self.toff = toff
		self.scale = scale
	
	def theoretical_max(self):
		"""
		Returns the absolute highest value of reproduction possible, achieved when repair is zero.
		This value is just equal to scale (as said before).

		Returns:
		float: The theoretical maximum fecundity.
		"""
		return self.scale
	
	def ret_val(self):
		"""
		Returns the peak reproductive value given the current repair investment
		
		Returns:
		float: The peak reproductive value
		"""
		return self.toff.ret_val() * self.scale


class PeakReproductiveValueMulti:
	"""
	An object that represents the highest reproductive value possible given the current investment in repair. This
	object can take multiple repair investments with different scales (such that the final value returned
	is scale_1*toff_1(r_1) + scale_2*toff_2(r_2) + ...)
	"""
	def __init__(self, toffs, scales):
		"""
		Init function

		Parameters:
		- toffs (list): List of trade off function objects (like ReproductiveScaling)
		- scales (list): List of the scales of reproduction, such that the absolute highest value of reproduction
						 theoretically possible with no repair is equal to the scale variable.
		"""
		self.toffs = toffs
		self.scales = scales
	
	def theoretical_max(self):
		"""
		Returns the absolute highest value of reproduction possible, achieved when all repair investments
		are equal to zero.

		Returns:
		float: The theoretical maximum fecundity.
		"""
		return sum(self.scales)
	
	def ret_val(self):
		"""
		Returns the peak reproductive value given the current repair investment

		Returns:
		float: The peak reproductive value
		"""
		cur_val = 0
		for i in range(len(self.toffs)):
			cur_val += self.toffs[i].ret_val() * self.scales[i]
		return cur_val


class MatAge:
	"""
	An object that represents the age of reproductive maturation and the first age of reproduction. This age is
	calculated from some min_age parameter and the peak_reproductive_value object. The equation used is
	min_age/(peak_reproductive_value/theoretical_max). For the default form of (1-r) for peak_reproductive_value, the
	equation reduces to min_age/(1-r)
	"""
	def __init__(self, min_age, peak_reproductive_value):
		"""
		Init function

		Parameters:
		- min_age (float): The minimum possible age of reproductive maturity
		- peak_reproductive_value (object): The peak reproductive value object (either PeakReproductiveValue or
			PeakReproductiveValueMulti)
		"""
		self.min_age = min_age
		self.peak_reproductive_value = peak_reproductive_value
	
	def ret_val(self):
		"""
		Returns the maturation age given the current repair investment

		Returns:
		float: The maturation age
		"""
		return self.min_age / (self.peak_reproductive_value.ret_val() / self.peak_reproductive_value.theoretical_max())


class ConstMatAge:
	"""
	An object that represents the age of reproductive maturation and the first age of reproduction. In this object,
	the maturation age is constant and does not depend on repair.
	"""
	def __init__(self, min_age):
		"""
		Init function

		Parameters:
		- min_age (float): The age of reproductive maturity
		"""
		self.min_age = min_age
	
	def ret_val(self):
		"""
		Returns the maturation age

		Returns:
		float: The maturation age
		"""
		return self.min_age


class DecayingRepro:
	"""
	An object that represents the instantaneous reproduction function, in which reproduction tends to decline over time
	at a rate determined by the rate of senescence.
	"""
	def __init__(self, decay_funcs, mat_age, peak_reproductive_value):
		"""
		Init function

		Parameters:
		- decay_funcs (list): List of monotonically decreasing function objects (like GompertzRepAMR or Extrinsic)
		- mat_age (MatAge or ConstMatAge): Maturation age object (defines when reproduction should start)
		- peak_reproductive_value (PeakReproductiveValue or PeakReproductiveValueMulti) peak reproductive value object
		"""
		self.decay_funcs = decay_funcs
		self.mat_age = mat_age
		self.peak_reproductive_value = peak_reproductive_value
	
	def ret_val_at_time(self, x):
		"""
		Calculates the instantaneous reproductive output

		Parameters:
		- x (float): The current age

		Returns:
		float: The current instantaneous fecundity
		"""
		if x < self.mat_age.ret_val():
			return 0
		else:
			cur_rep_level = self.peak_reproductive_value.ret_val()
			for func in self.decay_funcs:
				cur_rep_level *= func.integral_form(x, self.mat_age.ret_val(), decay_start=self.mat_age.ret_val())
			return cur_rep_level


class ConstRepro:
	"""
	An object that represents the instantaneous reproduction function. There is no reproductive senescence here.
	"""
	def __init__(self, mat_age, peak_reproductive_value):
		"""
		Init function

		Parameters:
		- decay_funcs (list): List of monotonically decreasing function objects (like GompertzRepAMR or Extrinsic)
		- mat_age (MatAge or ConstMatAge): Maturation age object (defines when reproduction should start)
		- peak_reproductive_value (PeakReproductiveValue or PeakReproductiveValueMulti) peak reproductive value object
		"""
		self.mat_age = mat_age
		self.peak_reproductive_value = peak_reproductive_value
	
	def ret_val_at_time(self, x):
		"""
		Calculates the instantaneous reproductive output (constant here past maturity)

		Parameters:
		- x (float): The current age

		Returns:
		float: The current instantaneous fecundity (0 if before maturation age, peak reproductive value otherwise)
		"""
		if x < self.mat_age.ret_val():
			return 0
		else:
			return self.peak_reproductive_value.ret_val()


class DecayingReproSizeGrowth:
	"""
	An object that represents the instantaneous reproduction function, in which reproduction tends to decline over time
	at a rate determined by the rate of senescence. However, reproduction also increases linearly as animals age,
	representing indeterminate growth.
	"""
	def __init__(self, decay_funcs, mat_age, peak_reproductive_value, ind_growth):
		"""
		Init function

		Parameters:
		- decay_funcs (list): List of monotonically decreasing function objects (like GompertzRepAMR or Extrinsic)
		- mat_age (MatAge or ConstMatAge): Maturation age object (defines when reproduction should start)
		- peak_reproductive_value (PeakReproductiveValue or PeakReproductiveValueMulti) peak reproductive value object
		- ind_growth (float): The rate of growth (multiplied by adult age and added to reproduction)
		"""
		self.decay_funcs = decay_funcs
		self.mat_age = mat_age
		self.peak_reproductive_value = peak_reproductive_value
		self.ind_growth = ind_growth
	
	def ret_val_at_time(self, x):
		"""
		Calculates the instantaneous reproductive output. Affected by both senescence and size growth.

		Parameters:
		- x (float): The current age

		Returns:
		float: The current instantaneous fecundity
		"""
		if x < self.mat_age.ret_val():
			return 0
		else:
			cur_rep_level = self.peak_reproductive_value.ret_val()
			for func in self.decay_funcs:
				cur_rep_level *= func.integral_form(x, self.mat_age.ret_val(), decay_start=self.mat_age.ret_val())
			time_since_mat = x-self.mat_age.ret_val()
			rep_add = time_since_mat * self.ind_growth
			cur_rep_level += rep_add
			return cur_rep_level


class DecayingReproSizeGrowthMul:
	"""
	An object that represents the instantaneous reproduction function, in which reproduction tends to decline over time
	at a rate determined by the rate of senescence. However, reproduction also increases as animals age,
	representing indeterminate growth. The ind_growth term is multiplied instead of added to reproduction.
	"""
	def __init__(self, decay_funcs, mat_age, peak_reproductive_value, ind_growth):
		"""
		Init function

		Parameters:
		- decay_funcs (list): List of monotonically decreasing function objects (like GompertzRepAMR or Extrinsic)
		- mat_age (MatAge or ConstMatAge): Maturation age object (defines when reproduction should start)
		- peak_reproductive_value (PeakReproductiveValue or PeakReproductiveValueMulti) peak reproductive value object
		- ind_growth (float): The rate of growth (multiplied by adult age and scaled to reproduction)
		"""
		self.decay_funcs = decay_funcs
		self.mat_age = mat_age
		self.peak_reproductive_value = peak_reproductive_value
		self.ind_growth = ind_growth
	
	def ret_val_at_time(self, x):
		"""
		Calculates the instantaneous reproductive output. Affected by both senescence and size growth.

		Parameters:
		- x (float): The current age

		Returns:
		float: The current instantaneous fecundity
		"""
		if x < self.mat_age.ret_val():
			return 0
		else:
			cur_rep_level = self.peak_reproductive_value.ret_val()
			for func in self.decay_funcs:
				cur_rep_level *= func.integral_form(x, self.mat_age.ret_val(), decay_start=self.mat_age.ret_val())
			time_since_mat = x - self.mat_age.ret_val()
			#rep_mul = 1+(time_since_mat * self.ind_growth)
			if time_since_mat > 50:
				time_since_mat = 50
			rep_mul = math.exp(time_since_mat*self.ind_growth)
			cur_rep_level *= rep_mul
			return cur_rep_level


class SurvFunc:
	"""
	An object that represents the survival function, l(x).
	"""
	def __init__(self, mort_funcs):
		"""
		Init function

		Parameters:
		- mort_funcs (list): List of monotonically decreasing function objects (like GompertzRepAMR or Extrinsic)
		"""
		self.mort_funcs = mort_funcs
	
	def ret_val_at_time(self, x):
		"""
		Calculates the fraction of individuals surviving to age x.

		Parameters:
		- x (float): The age

		Returns:
		float: The fraction of individuals (from 0 to 1) surviving to age x
		"""
		cur_surv = 1
		for func in self.mort_funcs:
			cur_surv *= func.integral_form(x, 0)
		return cur_surv
	
	def max_age(self):
		"""
		Finds an age at which a fraction of individuals less than 1e-10 survives. Returns 10000 if no age before
		this is found. Faster than age_at_frac_remaining(1e-10)

		Returns:
		float: The "max age" of the population
		"""
		cur_val = 50
		while self.ret_val_at_time(cur_val) > 1e-10 and cur_val < 10000:
			cur_val*=2
		return cur_val

	
	def age_at_frac_remaining(self, frac):
		"""
		Calculates the age at which `frac` individuals are left. For instance, age_at_frac_remaining(0.5) will
		give the median lifespan of the population.

		Parameters:
		- frac (float): The fraction of individuals (1 to 0) left surviving

		Returns:
		float: The age at which the population reaches `frac` individuals surviving
		"""
		obj_func_abs = lambda x: self.ret_val_at_time(x) - frac
		res = minimize(obj_func_abs, np.asarray([50]),bounds=[[0,None]], method='nelder-mead',
					   options={"fatol": 1E-4, "xatol": 1E-4})
		cur_val = res['x'][0]
		return cur_val


class SurvFuncDK:
	"""
	An object that represents the survival function, l(x). Uses the original Drenos and Kirkwood model that
	incorporates juvenile mortality as a separate parameter.
	"""
	def __init__(self, mort_funcs, mat_age, juv_surv):
		"""
		Init function

		Parameters:
		- mort_funcs (list): List of monotonically decreasing function objects (like GompertzRepAMR or Extrinsic)
		- mat_age (MatAge or ConstMatAge): Maturation age object
		- juv_surv (float): The fraction of newborns surviving to adulthood
		"""
		self.mort_funcs = mort_funcs
		self.mat_age = mat_age
		self.juv_surv = juv_surv
	
	def ret_val_at_time(self, x):
		"""
		Calculates the fraction of individuals surviving to age x.

		Parameters:
		- x (float): The age

		Returns:
		float: The fraction of individuals (from 0 to 1) surviving to age x
		"""
		a = self.mat_age.ret_val()
		cur_surv = self.juv_surv
		if x < a:
			return cur_surv
		for func in self.mort_funcs:
			cur_surv *= func.integral_form(x, a, decay_start=a)
		return cur_surv
	
	def max_age(self):
		"""
		Finds an age at which a fraction of individuals less than 1e-10 survives. Returns 10000 if no age before
		this is found. Faster than age_at_frac_remaining(1e-10)

		Returns:
		float: The "max age" of the population
		"""
		cur_val = 50
		while self.ret_val_at_time(cur_val) > 1e-10 and cur_val < 10000:
			cur_val*=2
		return cur_val

	def age_at_frac_remaining(self, frac):
		"""
		Calculates the age at which `frac` individuals are left. For instance, age_at_frac_remaining(0.5) will
		give the median lifespan of the population.

		Parameters:
		- frac (float): The fraction of individuals (1 to 0) left surviving

		Returns:
		float: The age at which the population reaches `frac` individuals surviving
		"""
		obj_func_abs = lambda x: self.ret_val_at_time(x) - frac
		res = minimize(obj_func_abs,np.asarray([50]),method='nelder-mead',options={"fatol":1E-8,"xatol":1E-8})
		cur_val = res['x'][0]
		return cur_val


class RSelPop:
	"""
	An object that represents a population under r-selection.
	"""
	def __init__(self, rep_func, surv_func, use_romberg=False):
		"""
		Init function

		Parameters:
		- rep_func (object): Reproductive function object (like DecayingRepro)
		- surv_func (object): Survival function object (like SurvFunc)
		- use_romberg (bool): Flag to use romberg package for integrating the Euler-Lotka equation (default False)
		"""
		self.rep_func = rep_func
		self.surv_func = surv_func
		self.use_romberg = use_romberg
	
	def ret_val(self):
		"""
		Calculates the rate of increase from the Euler-Lotka equation here (we call it g to avoid confusion with repair).

		Returns:
		float: The intrinsic rate of population growth `g` (known as r in other publications)
		"""
		max_age = self.surv_func.max_age()
		obj_func = lambda g: integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x)
												 * math.exp(-x * g), 0, max_age, use_romberg=self.use_romberg) - 1
		obj_func_abs = lambda x : abs(obj_func(x))
		if obj_func(0) < 0:
			return 0
		res = minimize(obj_func_abs,np.asarray([0.15]),bounds=[[0,None]],method='nelder-mead',options={"fatol":1E-6,"xatol":1E-6})
		if res['x'][0] == 0:
			res = minimize(obj_func_abs,np.asarray([0.001]),bounds=[[0,None]],method='nelder-mead',options={"fatol":1E-6,"xatol":1E-6})
		g_val = res['x'][0]
		if g_val == 0.15:
			return 0
		return(g_val)

	def hamilton_mort_additive(self,age,this_growth=None):
		"""
		The hamilton indicator for a mutation with an additive effect on mortality at a certain age.
		
		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for additive mortality
		"""
		max_age = self.surv_func.max_age()
		if age > max_age:
			return 0
		if this_growth is None:
			this_growth = self.ret_val()
		res = integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x) *
							math.exp(-this_growth * x), age,max_age) / integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x) *
							math.exp(-this_growth * x) * x, 0,max_age)
		return res
	
	def hamilton_mort_proportional(self,age,this_growth=None):
		"""
		The hamilton indicator for a mutation with a proportional effect on mortality at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for proportional mortality
		"""
		mort_val = 0
		for mort_func in self.surv_func.mort_funcs:
			mort_val += mort_func.ret_val_at_time(age)
		return self.hamilton_mort_additive(age,this_growth)*mort_val
	
	def hamilton_fert_additive(self,age,this_growth=None):
		"""
		The hamilton indicator for a mutation with an additive effect on fertility at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for additive fertility
		"""
		max_age = self.surv_func.max_age()
		if age > max_age:
			return 0
		if this_growth is None:
			this_growth = self.ret_val()
		numerator = math.exp(-this_growth * age) * self.surv_func.ret_val_at_time(age)
		res = integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x) * x *
								  math.exp(-this_growth*x), 0, max_age)
		return numerator / res
	
	def hamilton_fert_proportional(self, age,this_growth=None):
		"""
		The hamilton indicator for a mutation with a proportional effect on fertility at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for proportional fertility
		"""
		return self.rep_func.ret_val_at_time(age) * self.hamilton_fert_additive(age,this_growth)


class KSelPop:
	"""
	An object that represents a population under K-selection.
	"""
	def __init__(self, rep_func, surv_func, use_romberg=False):
		"""
		Init function

		Parameters:
		- rep_func (object): Reproductive function object (like DecayingRepro)
		- surv_func (object): Survival function object (like SurvFunc)
		- use_romberg (bool): Flag to use romberg package for integrating the Euler-Lotka equation (default False)
		"""
		self.rep_func = rep_func
		self.surv_func = surv_func
		self.use_romberg = use_romberg
	
	def ret_val(self):
		"""
		Calculates the lifetime reproductive output (R0), equal to the integral of l(x)m(x)

		Returns:
		float: R0
		"""
		max_age = self.surv_func.max_age()
		first_age = self.rep_func.mat_age.ret_val()
		cur_val = integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x),
							first_age, max_age, self.use_romberg)
		return cur_val
	
	def hamilton_mort_additive(self, age):
		"""
		The hamilton indicator for a mutation with an additive effect on mortality at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for additive mortality
		"""
		max_age = self.surv_func.max_age()
		if age > max_age:
			return 0
		res = integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x), age, max_age) / \
			  integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x) * x, 0, max_age)
		return res
	
	def hamilton_mort_proportional(self, age):
		"""
		The hamilton indicator for a mutation with a proportional effect on mortality at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for proportional mortality
		"""
		mort_val = 0
		for mort_func in self.surv_func.mort_funcs:
			mort_val += mort_func.ret_val_at_time(age)
		return self.hamilton_mort_additive(age) * mort_val
	
	def hamilton_fert_additive(self, age):
		"""
		The hamilton indicator for a mutation with an additive effect on fertility at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for additive fertility
		"""
		max_age = self.surv_func.max_age()
		if age > max_age:
			return 0
		numerator = self.surv_func.ret_val_at_time(age)
		res = integrate(lambda x: self.surv_func.ret_val_at_time(x) * self.rep_func.ret_val_at_time(x) * x, 0, max_age)
		return numerator / res
	
	def hamilton_fert_proportional(self, age):
		"""
		The hamilton indicator for a mutation with a proportional effect on fertility at a certain age.

		From Baudisch 2005 "Hamilton's indicators of the force of selection"

		Returns:
		float: Hamilton indicator for proportional fertility
		"""
		return self.rep_func.ret_val_at_time(age) * self.hamilton_fert_additive(age)


def eval_growth(growth_func, rep_vars, inputs):
	"""
	Helper function that returns the population growth (either g or R0) by assigning the repair values to the
	repair variables

	Parameters:
	- growth_func (object): The population object (like RSelPop or KSelPop)
	- rep_vars (list): List of Repair objects
	- inputs (list): List of float values to be assigned to rep_vars

	Returns:
	float: The ret_val of the growth_func when rep_vars are set to inputs
	"""
	for i,val in enumerate(inputs):
		rep_vars[i].val = val
	return growth_func.ret_val()


def optimize_growth(growth_func, repair_vars, repair_ranges, start_point=None):
	"""
	Tries to optimize the repair investments to maximize population growth (or R0)

	Parameters:
	- growth_func (object): The population object (like RSelPop or KSelPop)
	- rep_vars (list): List of Repair objects
	- repair_ranges (list): List of lists dictating the lower and upper bound of the investment search space
							for each repair variable (e.g. [[0,0.8],[0,0.5]]

	Returns:
	float: max_growth, max_point: Maximum value of g or R0, and the maximum values of the repair variables which
									achieve it.
	"""
	if start_point is None:
		max_growth = -float("inf")
		max_point = None
		point_list = itertools.product(*[np.linspace(x[0] + 0.01, x[1] - 0.0001, num=30) for x in repair_ranges])
		for points in point_list:
			for i, point in enumerate(points):
				repair_vars[i].val = point
			this_growth = growth_func.ret_val()
			if this_growth > max_growth:
				max_growth = this_growth
				max_point = points
		start_point = max_point
	
	max_point = start_point
	max_growth = eval_growth(growth_func,repair_vars,max_point)
	res = minimize(lambda x: -eval_growth(growth_func, repair_vars, x),
				   np.asarray(start_point), bounds=repair_ranges,
				   method='nelder-mead',
				   options={"xatol":1e-5,"fatol":1e-4})
	if -res['fun'] > max_growth:
		max_growth = -res['fun']
		max_point = list(res['x'])
	if max_growth == 0:
		res = scipy.optimize.brute(lambda x: -eval_growth(growth_func, repair_vars, x),
					   repair_ranges, Ns=20, full_output=True,finish=None)
		max_growth = -res[1]
		max_point = res[0]
		res = minimize(lambda x : -eval_growth(growth_func,repair_vars,x), np.asarray(max_point), bounds=repair_ranges,method='nelder-mead',options={})
		max_growth = -res['fun']
		max_point = list(res['x'])
	for i in range(len(max_point)):
		new_point = [x for x in max_point]
		new_point[i] = repair_ranges[i][1]
		new_growth = eval_growth(growth_func,repair_vars,new_point)
		if new_growth > max_growth:
			max_growth = new_growth
			max_point = new_point
	fin_test = [x[1] for x in repair_ranges]
	fin_growth = eval_growth(growth_func,repair_vars,fin_test)
	if fin_growth > max_growth:
		max_growth = fin_growth
		max_point = fin_test
	max_growth = eval_growth(growth_func, repair_vars, max_point)
	return (max_growth,max_point)
	
