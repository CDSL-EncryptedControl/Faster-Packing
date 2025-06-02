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
		LogP:    []int{51},
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

	// /// ============== Plant model ==============
	// A := [][]float64{
	// 	{1.64229732608657, 3.53029220299480, 2.25830940585494, -1.66913412899748, -0.826587098687241},
	// 	{1.38573096400922, -0.690427306343379, -0.276045798482360, 0.526163274322444, 0.160978063182478},
	// 	{0.844228121278874, 0.195353485340086, -0.0790399929447087, -0.113711975486802, -0.0830905241387200},
	// 	{12.8867357478498, 9.29152468776712, 5.92466342682760, -6.31980680583488, -2.24050901566593},
	// 	{-0.768103327549902, -0.260198593549999, -0.0212707258138640, 0.173241959148395, 0.107190092415570},
	// }
	// B := [][]float64{
	// 	{12.3405535089277, 0.0384208786823148},
	// 	{0.00691726599398608, 0.00615342592981413},
	// 	{-15.9998195696488, 0.000489666760710950},
	// 	{0.00376477989315791, -0.507076963222898},
	// 	{0.000284391490834563, -0.952939300459448},
	// }
	// C := [][]float64{
	// 	{2.51076385209752, 1.83468871140524, 1.18288134573927, -1.23701748848539, -0.639935268982959},
	// 	{-1.09907974591350, -0.378353541288560, -0.0346764721467886, 0.242599379784725, 0.150575344095704},
	// }
	// // Plant initial state
	// xp_ini := []float64{
	// 	0,
	// 	0,
	// 	0.1,
	// 	-0.1,
	// 	0.6,
	// }

	// // ============== Pre-designed controller ==============
	// // F must be an integer matrix
	// F := [][]float64{
	// 	{3, 1, 0, 0, 0},
	// 	{0, 0, 1, 0, 0},
	// 	{-4, 0, 0, 0, 0},
	// 	{0, 0, 0, 1, 1},
	// 	{0, 0, 0, 2, 0},
	// }
	// G := [][]float64{
	// 	{2, 3},
	// 	{-1.2, -4},
	// 	{-0.1, -1},
	// 	{5, -0.3},
	// 	{0, 0.7},
	// }
	// H := [][]float64{
	// 	{0.25, 0, 0, 0, 0},
	// 	{0, 0, 0, -2.1, 0},
	// }
	// // Controller initial state
	// x_ini := []float64{
	// 	0,
	// 	0,
	// 	0,
	// 	0,
	// 	0,
	// }

	/// ============== Plant model ==============
	/// ============== Plant model ==============
	A := [][]float64{
	{-4.9535, -1.3701, 2.0157, 1.0929},
	{5.4838, 3.0300, -3.8440, -1.9888},
	{0.9319, 0.5722, -1.2467, 0.5866},
	{-2.4378, -0.9447, -2.0371, -0.6299},
	}	
	B := [][]float64{
	{1.3993, -0.0344},
	{2.0586, -0.0405},
	{0.0968, 0.4669},
	{0.1186, 0.6871},
	}
	C := [][]float64{
	{-0.5224, -0.2219, -0.3423, -0.1006},
	{-0.9765, -0.5500, 0.8802, 0.4234},
	}
	// Plant initial state
	xp_ini := []float64{
		0,
		0,
		0.1,
		-0.1,
	}

	// ============== Pre-designed controller ==============
	// F must be an integer matrix
	F := [][]float64{
	{1, 1, 0, 0},
	{2, 0, 0, 0},
	{0, 0, 1, 1},
	{0, 0, 2, 0},
	}
	
	G := [][]float64{
	{2.7, 3.2},
	{-1.3, -4.9},
	{-0.1000, -1.0000},
	{5.0000, -0.3000},
	}
	
	H := [][]float64{
		{1, 0, 0, 0},
		{0, 0, 3, 0},
	}
	// Controller initial state
	x_ini := []float64{
		0,
		0,
		0,
		0,
	}


	// dimensions
	n := len(F)
	m := len(H)
	p := len(G[0])

	// ============== Quantization parameters ==============
	s := 1 / 10000.0  // 1e+4
	L := 1 / 10000000000.0 // 1e+10
	r := 1 / 100000.0 // 1e+5
	fmt.Printf("Scaling parameters 1/L: %v, 1/s: %v, 1/r: %v \n", 1/L, 1/s, 1/r)
	// *****************************************************************
	// *****************************************************************

	// ============== Encryption settings ==============
	// Set parameters
	levelQ := params.QCount() - 1
	levelP := params.PCount() - 1
	ringQ := params.RingQ()

	// Compute tau
	// least power of two greater than n, p_, and m
	maxDim := math.Max(math.Max(float64(n), float64(m)), float64(p))
	tau := int(math.Pow(2, math.Ceil(math.Log2(maxDim))))

	// Generate monomials for unpack
	logn := int(math.Log2(float64(tau)))
	monomials := make([]ring.Poly, logn)
	for i := 0; i < logn; i++ {
		monomials[i] = ringQ.NewPoly()
		idx := params.N() - params.N()/(1<<(i+1))
		monomials[i].Coeffs[0][idx] = 1
		ringQ.MForm(monomials[i], monomials[i])
		ringQ.NTT(monomials[i], monomials[i])
	}

	// tau = 8
	// galEls = 3개짜리
	// galEls = [8+1, 4+1, 2+1]

	// Generate Galois elements for unpack
	galEls := make([]uint64, int(math.Log2(float64(tau))))
	for i := 0; i < int(math.Log2(float64(tau))); i++ {
		galEls[i] = uint64(tau/int(math.Pow(2, float64(i))) + 1)
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
	// HBar := utils.ScalMatMult(1/s, H)
	HBar := utils.ScalMatMult(1, H)

	// Encryption
	// Dimension: 1-by-(# of columns)
	ctF := RGSW.EncPack(F, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctG := RGSW.EncPack(GBar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)
	ctH := RGSW.EncPack(HBar, tau, encryptorRGSW, levelQ, levelP, ringQ, params)

	// ============== Simulation ==============
	// Number of simulation steps
	iter := 200
	fmt.Printf("Number of iterations: %v\n", iter)

	// *****************
	// 1) Plant + unencrypted (original) controller
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
	xpEnc = append(xpEnc, xp_ini)

	// Plant state
	xp = xp_ini

	// Dimension: 1-by-(# of elements)
	xBar := utils.ScalVecMult(1/(r*s), x_ini)
	xCtPack := RLWE.EncPack(xBar, tau, 1/L, *encryptorRLWE, ringQ, params)

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
		yCtPack := RLWE.EncPack(yBar, tau, 1/L, *encryptorRLWE, ringQ, params)

		// **** Encrypted Controller ****
		// Unpack state
		xCt := RLWE.UnpackCt(xCtPack, n, tau, evaluatorRLWE, ringQ, monomials, params)

		// Unpack input
		yCt := RLWE.UnpackCt(yCtPack, p, tau, evaluatorRLWE, ringQ, monomials, params)

		// Compute output
		uCtPack := RGSW.MultPack(xCt, ctH, evaluatorRGSW, ringQ, params)

		// **** Actuator ****
		// Decrypt and Unapck
		// u := RLWE.DecUnpack(uCtPack, m, tau, *decryptorRLWE, r*s*s*L, ringQ, params)
		u := RLWE.DecUnpack(uCtPack, m, tau, *decryptorRLWE, s*L, ringQ, params)

		// **** Encrypted Controller ****
		// State update
		FxCt := RGSW.MultPack(xCt, ctF, evaluatorRGSW, ringQ, params)
		GyCt := RGSW.MultPack(yCt, ctG, evaluatorRGSW, ringQ, params)
		xCtPack = RLWE.Add(FxCt, GyCt, params)

		period[i] = []float64{float64(time.Since(startPeriod[i]).Microseconds()) / 1000}

		// **** Plant ****
		// State update
		xp = utils.VecAdd(utils.MatVecMult(A, xp), utils.MatVecMult(B, u))

		// Save data
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
	utils.DataExport(uDiff, "./uDiff_MIMO_13.csv")
	
	// Elapsed time
	utils.DataExport(period, "./period.csv")
	utils.DataExport(period, "./period_MIMO_13.csv")

	// nominal trajectory
	utils.DataExport(xcUnenc, "./nomctrlstate.csv")
	utils.DataExport(xpUnenc, "./nomstate.csv")

	// nominal input
	utils.DataExport(uUnenc, "./nominput.csv")
}
