# CDSL: Cryptography for Dynamic Systems Library


CDSL provides codes for implementing secure dynamic systems based on modern cryptography.
The library features linear dynamic controllers operating over homomorphically encrypted data implemented using [Lattigo](https://github.com/tuneinsight/lattigo) version 6.1.0.
The encrypted controllers are designed based on the state-of-the-art methods developed by CDSL, [SNU](https://post.cdsl.kr/) and [SEOULTECH](https://junsookim4.wordpress.com/).


---

### Features


Given a plant 

$$
\begin{aligned}
x_p(t+1) &= Ax_p(t) + Bu(t), \quad x_p(0) = x_p^{\mathsf{ini}} \\
y(t) &= Cx_p(t)
\end{aligned}
$$

and a pre-designed stabilizng controller (which is controllable and observable)

$$
\begin{aligned}
x(t+1) &= Fx(t) + Gy(t), \quad x(0) = x^{\mathsf{ini}} \\
u(t) &= Hx(t)
\end{aligned}
$$

this code provides two methods to operate the linear dynamic controller over encrypted data, using a Ring-LWE based cryptosystem. 


- `ctrRGSW` [1]: Supports unlimited number of recursive homomorhpic multiplications without the use of bootstrapping. More specifically, the encrypted controller state is recursively multiplied to the encrypted state matrix without decryption. The effect of error growth is suppressed by the stability of the closed-loop system. 
    - `ctrRGSW/noPacking`: Naive implementation that does not use packing. 
    - `ctrRGSW/packing`: A novel "coefficient packing" technique is applied, resulting in enhanced computation speed and memory efficiency   
    - `ctrRGSW/conversion.m`: Converts the state matrix of the controller into integers based on the apporach of [3]:
       - Given $F$ and $H$, it finds an appropriate $R$ such that $F-RH$ is an integer matrix. Then, the state dynamics of the controller can be rewritten as follows regarding $u(t)$ as a fed-back input.
       - One may use other methods, such as [4-6], to convert the state matrix into integers while preserving the control performance without using re-encryption. 


$$
x(t+1) = (F-RH)x(t) + Gy(t) + Ru(t)
$$




- `ctrRLWE` [2]: 

- In this example, we use the four-tank system [7] as the plant and pre-design an observer-based controller. Please refer to [8] for more details regarding the setting.

---

### How to use
1. Download or clone this repository using
```
git clone https://github.com/CDSL-EncryptedControl/CDSL.git
```


2. Change the directory to a folder you wish to use and run the `main.go` file. For example,

```
cd ctrRGSW/noPacking
go run main.go  
```
or
```
cd ctrRGSW/packing
go run main.go  
```
or
```
cd ctrRLWE
go run main.go  
```

---

### License
Distributed under the MIT License. See `LICENSE` for more information.

---

### Contributing

We **greatly appreciate** any typo and factual corrections, no matter how small the change is.

If you spot an error in our code, please report it to us without hesitation.
If the issue is not completely obvious, please provide justifications and details as necessary. 

---


### Contact

Yeongjun Jang - jangyj@cdsl.kr

Joowon Lee - jwlee@cdsl.kr

Junsoo Kim - junsookim@seoultech.ac.kr

---

### Acknowledgements
- This work was supported by the National Research Foundation of Korea(NRF) grant funded by the Korea government(MSIT) (No. RS-2024-00353032).
- Special thanks to [Seonhong Min](https://snu-lukemin.github.io/), [Hyesun Kwak](https://hyesunkwak.github.io/), and [Yongsoo Song](https://yongsoosong.github.io/) with the Department of Computer Science and Engineering, Seoul National University, for the great help.

---

### References
[1] [Y. Jang, J. Lee, S. Min, H. Kwak, J. Kim, and Y. Song, "Ring-LWE based encrypted controller with unlimited number of recursive multiplications and effect of error growth," 2024, arXiv:2406.14372.](https://arxiv.org/abs/2406.14372)

[2] [J. Lee, D. Lee, J. Kim, and H. Shim, "Encrypted dynamic control exploiting limited number of multiplications and a method using RLWE-based cryptosystem," _IEEE Trans. Syst. Man. Cybern.: Syst._, vol. 55, no. 1, pp. 158-169, 2025.](https://ieeexplore.ieee.org/abstract/document/10730788)

[3] [J. Kim, H. Shim, and K. Han, "Dynamic controller that operates over homomorphically encrypted data for infinite time horizon," _IEEE Trans. Autom. Control_, vol. 68, no. 2, pp. 660-672, 2023.](https://ieeexplore.ieee.org/abstract/document/9678042)

[4] [J. Kim, H. Shim, H. Sandberg, and K. H. Johansson, “Method for running dynamic systems over encrypted data for infinite time horizon without bootstrapping and re-encryption,” in Proc. 60th IEEE Conf. Decision Control, 2021, pp. 5614–5619.](https://ieeexplore.ieee.org/abstract/document/9682828?casa_token=LHR79rToQ7oAAAAA:Wz1AzFWR7VW6DYKUhLFYcoXtpMx4AIT9E_krpOpFy7QUO5lSkvPf_0ZZgPsdp65ZzaGx-ejlPA)

[5] [M. S. Tavazoei, “Non-minimality of the realizations and possessing state matrices with integer elements in linear discrete-time controllers,” IEEE Trans. Autom. Control, vol. 68, no. 6, pp. 3698–3703, 2023.](https://ieeexplore.ieee.org/abstract/document/9835020?casa_token=_rdGjQLc7ZEAAAAA:QLxzC1QlnNVYriMTL1gbSjtv5U2oTwfVO5OqVFfGS0Qpz8hx7exSuJKJ9H8XBh_qDucoZt8oBg)

[6] [J. Lee, D. Lee, S. Lee, J. Kim, and H. Shim, “Conversion of controllers to have integer state matrix for encrypted control: Non-minimal order approach,” in Proc. 62nd IEEE Conf. Decision Control, 2023, pp. 5091–5096.](https://ieeexplore.ieee.org/abstract/document/10383200?casa_token=lbob37tAZ-MAAAAA:vAVUmuIngRzHYefqaYHQM5TfukcAI7Lh1YmYngqcLYMj74Mtzq0xGybkntfWSd-DKwogxrvnxg)

[7] [K. Johansson, “The quadruple-tank process: A multivariable laboratory process with an adjustable zero,” IEEE Trans. Control Sys. Technol., vol. 8, no. 3, pp. 456–465, 2000.](https://ieeexplore.ieee.org/abstract/document/845876?casa_token=1CWEIgmKIscAAAAA:Hh3D4_xn5B8MWVoMpQHof8glwtWpGXMuddehBoKXbZAOh2WwsDlemeiWeZ6nAwQGThjhYYw1wQ)

[8] [Y. Jang, J. Lee, and J. Kim, "Documentation on encrypted dynamic control simulation code using Ring-LWE based cryptosystems," ](link)
