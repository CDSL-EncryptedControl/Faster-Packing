package main

import (
	"fmt"
	"math"
	"time"

	utils "github.com/CDSL-EncryptedControl/CDSL/utils"
	RGSW "github.com/CDSL-EncryptedControl/CDSL/utils/core/RGSW"
	RLWE "github.com/CDSL-EncryptedControl/CDSL/utils/core/RLWE"
	"github.com/tuneinsight/lattigo/v6/core/rgsw"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func main() {
	// *****************************************************************
	// ************************* User's choice *************************
	// *****************************************************************
	// ============== Encryption parameters ==============
	// Refer to ``Homomorphic encryption standard''
	params, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{
		// log2 of polynomial degree
		LogN: 13,
		// Size of ciphertext modulus (Q)
		LogQ: []int{56},
		// Size of special modulus (P)
		LogP: []int{51},
		NTTFlag: true,
	})
	fmt.Println("Degree of polynomials:", params.N())
	fmt.Println("Ciphertext modulus:", params.QBigInt())
	fmt.Println("Special modulus:", params.PBigInt())
	// Default secret key distribution
	// Each coefficient in the polynomial is uniformly sampled in [-1, 0, 1]
	fmt.Println("Secret key distribution (Ternary):", params.Xs())
	// Default error distribution
	// Each coefficient in the polynomial is sampled according to a
	// discrete Gaussian distribution with standard deviation 3.2 and bound 19.2
	fmt.Println("Error distribution (Discrete Gaussian):", params.Xe())

	/// Choose an example
	
	example := "numeric" // "invPen" or "numeric"
	fmt.Println("Example: ", example)
	
	var A, B, C [][]float64
    var xp_ini, x_ini []float64
    var F, G, H [][]float64
    var kappa int
    var ri []int
	
	if example == "invPen" {// inverted pendulum example
		/// ============== Plant model ==============
		A = [][]float64{
			{1,	0.0497869485895651,	0.00240013399977883,	3.99192178318363e-05},
			{0,	0.991480316476010,	0.0965381935679802,	0.00240013399977883},
			{0,	-0.000612279081576232,	1.04210214304071,	0.0506998325801739},
			{0,	-0.0246270901959133,	1.69541872243910,	1.04210214304071},
		}
		B = [][]float64{
			{0.00213051410434883},
			{0.0851968352399007},
			{0.00612279081576232},
			{0.246270901959133},
		}
		C = [][]float64{
			{1, 0, 0, 0},
		}
		// Plant initial state
		xp_ini = []float64{
			0,
			0,
			0.1,
			-0.1,
		}

		// ============== Pre-designed controller ==============
		// F must be an integer matrix
		F = [][]float64{
			{1, 1, 0, 0, 0, 0, 0, 0},
			{13, 0, 1, 0, 0, 0, 0, 0},
			{4, 0, 0, 1, 0, 0, 0, 0},
			{-10, 0, 0, 0, 1, 0, 0, 0},
			{0, 0, 0, 0, 0, 1, 0, 0},
			{0, 0, 0, 0, 0, 0, 1, 0},
			{0, 0, 0, 0, 0, 0, 0, 1},
			{0, 0, 0, 0, 0, 0, 0, 0},
		}
		
		kappa = 1
		
		ri = []int{
			0,
			}
		
		G = [][]float64{
			{-640.464761022051},
			{1715.36017949169},
			{-1489.05705877506},
			{389.491892102595},
			{27.2282414240346},
			{-0.604658260972883},
			{-2.36399144578786},
			{0.478391053561677},
		}
		H = [][]float64{
			{10, 0, 0, 0, 0, 0, 0, 0},
		}
		// Controller initial state
		x_ini = []float64{
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
		}
	} else if example == "numeric" { // numerical example
		/// ============== Plant model ==============
		A = [][]float64{
		{-4.9535, -1.3701, 2.0157, 1.0929},
		{5.4838, 3.0300, -3.8440, -1.9888},
		{0.9319, 0.5722, -1.2467, 0.5866},
		{-2.4378, -0.9447, -2.0371, -0.6299},
		}	
		B = [][]float64{
		{1.3993, -0.0344},
		{2.0586, -0.0405},
		{0.0968, 0.4669},
		{0.1186, 0.6871},
		}
		C = [][]float64{
		{-0.5224, -0.2219, -0.3423, -0.1006},
		{-0.9765, -0.5500, 0.8802, 0.4234},
		}
		// Plant initial state
		xp_ini = []float64{
			0,
			0,
			0.1,
			-0.1,
		}

		// ============== Pre-designed controller ==============
		// F must be an integer matrix
		F = [][]float64{
		{1, 1, 0, 0},
		{2, 0, 0, 0},
		{0, 0, 1, 1},
		{0, 0, 2, 0},
		}

		kappa = 2
		
		ri = []int{
		0, 
		2, 
		}
		
		G = [][]float64{
		{2.7, 3.2},
		{-1.3, -4.9},
		{-0.1000, -1.0000},
		{5.0000, -0.3000},
		}
		H = [][]float64{
			{1, 0, 0, 0},
			{0, 0, 3, 0},
		}
		// Controller initial state
		x_ini = []float64{
			0,
			0,
			0,
			0,
		}
	} else { // Default case
		A = [][]float64{
		{0},
		}	
		B = [][]float64{
		{0},
		}
		C = [][]float64{
		{0},
		}
		// Plant initial state
		xp_ini = []float64{
			0,
		}

		// ============== Pre-designed controller ==============
		// F must be an integer matrix
		F = [][]float64{
		{1},
		}

		kappa = 1
		
		ri = []int{
		0, 
		}
		
		G = [][]float64{
		{0},
		}
		H = [][]float64{
			{0},
		}
		// Controller initial state
		x_ini = []float64{
			0,
		}
	}
	
	// dimensions
	n := len(F)
	m := len(H)	

	// ============== Quantization parameters ==============
	s := 1 / 10000.0  // 1e+4
	L := 1 / 10000000000.0 // 1e+10
	r := 1 / 100000.0 // 1e+5
	fmt.Printf("Resolution level 1/r: %v \n", 1/r)
	fmt.Printf("Scaling parameters 1/L: %v, 1/s: %v \n", 1/L, 1/s)
	// *****************************************************************
	// *****************************************************************

	// ============== Encryption settings ==============
	// Set parameters
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQ := params.RingQ()

	// Compute tau
	// least power of two greater than n, p_, and m
	tau := int(math.Pow(2, math.Ceil(math.Log2(float64(m)))))
	
	// Define monomials
	monomial := ringQ.NewPoly() // X^{-N/n}
	monomial.Coeffs[0][params.N()-params.N()/n] = params.Q()[0] - 1
	ringQ.MForm(monomial, monomial)
	ringQ.NTT(monomial, monomial)
	
	monomials := make([]ring.Poly, kappa)
	monomials[0] = ringQ.NewPoly()
	monomials[0].Coeffs[0][0] = 1
	ringQ.MForm(monomials[0], monomials[0])
	ringQ.NTT(monomials[0], monomials[0])
	
	for i := 1; i < kappa; i++ {
	monomials[i] = ringQ.NewPoly()
	monomials[i].Coeffs[0][params.N()-ri[i]*params.N()/n] = params.Q()[0] - 1
	ringQ.MForm(monomials[i], monomials[i])
	ringQ.NTT(monomials[i], monomials[i])
	}
	
	// Generate Galois elements for unpack
	// n+1 n/2 + 1 ... 2 + 1
	galEls := make([]uint64, int(math.Log2(float64(n*tau))))
	for i := 0; i < int(math.Log2(float64(n*tau))); i++ {
		galEls[i] = uint64((n*tau)/int(math.Pow(2, float64(i))) + 1)
	}

	// Generate keys
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evkRGSW := rlwe.NewMemEvaluationKeySet(rlk)
	evkRLWE := rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...)

	// Define encryptor and evaluator
	encryptorRLWE := rlwe.NewEncryptor(params, sk)
	decryptorRLWE := rlwe.NewDecryptor(params, sk)
	encryptorRGSW := rgsw.NewEncryptor(params, sk)
	evaluatorRGSW := rgsw.NewEvaluator(params, evkRGSW)
	evaluatorRLWE := rlwe.NewEvaluator(params, evkRLWE)

	// ==============  Encryption of controller ==============
	// Quantization
	GBar := utils.ScalMatMult(1/s, G)
	HBar := utils.ScalMatMult(1, H)

	// Encryption
	// Dimension: 1-by-(# of columns)
	ctF := make([]*rgsw.Ciphertext, kappa)
	for i := 0; i < kappa; i++ {
		ctF[i] = RGSW.EncFi(F, ri[i], encryptorRGSW, levelQ, levelP, ringQ, params)	
	}
	ctG := RGSW.EncG(GBar, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctH := RGSW.EncH(HBar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)

	// ============== Simulation ==============
	// Number of simulation steps
	iter := 200
	fmt.Printf("Number of iterations: %v\n", iter)

	// *****************
	// 1) Plant + unencrypted (quantized) controller
	// *****************

	// State and output storage
	yUnenc := [][]float64{}
	uUnenc := [][]float64{}
	xcUnenc := [][]float64{}
	xpUnenc := [][]float64{}

	xpUnenc = append(xpUnenc, xp_ini)
	xcUnenc = append(xcUnenc, x_ini)

	// Plant state
	xp := xp_ini
	// Controller state
	x := utils.ScalVecMult(1/(r*s), x_ini)

	for i := 0; i < iter; i++ {
		y := utils.MatVecMult(C, xp)
		yBar := utils.RoundVec(utils.ScalVecMult(1/r, y))
		yBar = utils.ScalVecMult(r, yBar)
		u := utils.MatVecMult(H, x)
		uBar := utils.ScalVecMult(s, u)
		xp = utils.VecAdd(utils.MatVecMult(A, xp), utils.MatVecMult(B, uBar))
		x = utils.VecAdd(utils.MatVecMult(F, x), utils.MatVecMult(GBar, yBar))

		yUnenc = append(yUnenc, y)
		uUnenc = append(uUnenc, uBar)
		xcUnenc = append(xcUnenc, utils.ScalVecMult(r*s, x))
		xpUnenc = append(xpUnenc, xp)
	}

	// *****************
	// 2) Plant + encrypted controller
	// *****************

	// State and output storage
	yEnc := [][]float64{}
	uEnc := [][]float64{}
	xpEnc := [][]float64{}
	xcEnc := [][]float64{} // for debug
	xpEnc = append(xpEnc, xp_ini)
	xcEnc = append(xcEnc, x_ini)
	// Plant state
	xp = xp_ini

	// Dimension: 1-by-(# of elements)
	xBar := utils.ScalVecMult(1/(r*s), x_ini)
	xCtPack := RLWE.EncPack(xBar, n, 1/L, *encryptorRLWE, ringQ, params)

	// For time check
	period := make([][]float64, iter)
	startPeriod := make([]time.Time, iter)

	for i := 0; i < iter; i++ {
		// **** Sensor ****
		// Plant output
		y := utils.MatVecMult(C, xp)

		startPeriod[i] = time.Now()

		// Quantize and encrypt
		yBar := utils.RoundVec(utils.ScalVecMult(1/r, y))
		yBar = utils.ScalVecMult(r, yBar)
		yCtPack := RLWE.Ency(yBar, 1/L, *encryptorRLWE, ringQ, params)

		// **** Encrypted Controller ****

		// Compute output
		xTraceOutput := RLWE.Trace(xCtPack, n*tau, n, evaluatorRLWE, ringQ, params)
		uCtPack := RGSW.Mult(xTraceOutput, ctH, evaluatorRGSW, ringQ, params)

		// **** Actuator ****
		// Decrypt and Unpack
		u := RLWE.DecUnpack(uCtPack, m, n*tau, *decryptorRLWE, s*L, ringQ, params)

		// **** Encrypted Controller ****
		// State update
		FxCt := rlwe.NewCiphertext(params, xCtPack.Degree(), xCtPack.Level())
		ringQ.MulCoeffsMontgomery(xCtPack.Value[0], monomial, FxCt.Value[0])
		ringQ.MulCoeffsMontgomery(xCtPack.Value[1], monomial, FxCt.Value[1])
		
		for ii := 0; ii < kappa; ii++{
			xctTmp := rlwe.NewCiphertext(params, xCtPack.Degree(), xCtPack.Level())
			// xctTmp1 := rlwe.NewCiphertext(params, xCtPack.Degree(), xCtPack.Level())
			xTraceState := rlwe.NewCiphertext(params, xCtPack.Degree(), xCtPack.Level())
		
			ringQ.MulCoeffsMontgomery(xCtPack.Value[0], monomials[ii], xctTmp.Value[0])
			ringQ.MulCoeffsMontgomery(xCtPack.Value[1], monomials[ii], xctTmp.Value[1])
			xTraceState = RLWE.Trace(xctTmp, n, 1, evaluatorRLWE, ringQ, params)
			xctTmp1 := RGSW.Mult(xTraceState, ctF[ii], evaluatorRGSW, ringQ, params)
			
			ringQ.Add(FxCt.Value[0], xctTmp1.Value[0], FxCt.Value[0])
			ringQ.Add(FxCt.Value[1], xctTmp1.Value[1], FxCt.Value[1])
		}
		
		GyCt := RGSW.Mult(yCtPack, ctG, evaluatorRGSW, ringQ, params)
		xCtPack = RLWE.Add(FxCt, GyCt, params)

		period[i] = []float64{float64(time.Since(startPeriod[i]).Microseconds()) / 1000}

		// **** Plant ****
		// State update
		xp = utils.VecAdd(utils.MatVecMult(A, xp), utils.MatVecMult(B, u))

		// Save data
		// xc := RLWE.DecUnpack(xCtPack, n, n, *decryptorRLWE, r*s*L, ringQ, params)
		// xcEnc = append(xcEnc, xc)
		yEnc = append(yEnc, y)
		uEnc = append(uEnc, u)
		xpEnc = append(xpEnc, xp)
	}
	avgPeriod := utils.Average(utils.MatToVec(period))
	fmt.Println("Average elapsed time for a control period:", avgPeriod, "ms")

	// Compare plant input between 1) and 2)
	uDiff := make([][]float64, iter)
	for i := range uDiff {
		uDiff[i] = []float64{utils.Vec2Norm(utils.VecSub(uUnenc[i], uEnc[i]))}
	}

	// Export data ===============================================================

	// =========== Export data ===========

	// Plant state equipped with encrypted controller
	utils.DataExport(xpEnc, "./state.csv")

	// Plant intput from encrypted controller
	utils.DataExport(uEnc, "./uEnc.csv")

	// Plant output with encrypted controller
	utils.DataExport(yEnc, "./yEnc.csv")

	// Performance of encrypted controller
	utils.DataExport(uDiff, "./uDiff.csv")
	
	// Elapsed time
	utils.DataExport(period, "./period.csv")

	// nominal trajectory
	utils.DataExport(xcUnenc, "./nomctrlstate.csv")
	utils.DataExport(xpUnenc, "./nomstate.csv")

	// nominal input
	utils.DataExport(uUnenc, "./nominput.csv")
}
