#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <alsa/asoundlib.h>
#include <wiringPi.h>
#include <string.h>
#include <netinet/in.h>
#include <sys/ioctl.h>
#include <net/if.h>

#define PCM_DEVICE "default"

/////////////////////////////////////

int cols[12] = {23, 22, 21, 29, 28, 27, 26, 1, 4, 5, 6, 7};
int rows[5] = {0, 2, 3, 25, 24};


/************************************************
* FFT code from the book Numerical Recipes in C *
* Visit www.nr.com for the licence.             *
************************************************/

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <math.h>

#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

/*
 FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)

 Inputs:
	data[] : array of complex* data points of size 2*NFFT+1.
		data[0] is unused,
		* the n'th complex number x(n), for 0 <= n <= length(x)-1, is stored as:
			data[2*n+1] = real(x(n))
			data[2*n+2] = imag(x(n))
		if length(Nx) < NFFT, the remainder of the array must be padded with zeros

	nn : FFT order NFFT. This MUST be a power of 2 and >= length(x).
	isign:  if set to 1, 
				computes the forward FFT
			if set to -1, 
				computes Inverse FFT - in this case the output values have
				to be manually normalized by multiplying with 1/NFFT.
 Outputs:
	data[] : The FFT or IFFT results are stored in data, overwriting the input.
*/

void four1(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}

/********************************************************
* The following is a test routine that generates a ramp *
* with 10 elements, finds their FFT, and then finds the *
* original sequence using inverse FFT                   *
********************************************************/

// int main(int argc, char * argv[])
// {
// 	int i;
// 	int Nx;
// 	int NFFT;
// 	double *x;
// 	double *X;

// 	/* generate a ramp with 10 numbers */
// 	Nx = 10;
// 	printf("Nx = %d\n", Nx);
// 	x = (double *) malloc(Nx * sizeof(double));
// 	for(i=0; i<Nx; i++)
// 	{
// 		x[i] = i;
// 	}

// 	/* calculate NFFT as the next higher power of 2 >= Nx */
// 	NFFT = (int)pow(2.0, ceil(log((double)Nx)/log(2.0)));
// 	printf("NFFT = %d\n", NFFT);

// 	/* allocate memory for NFFT complex numbers (note the +1) */
// 	X = (double *) malloc((2*NFFT+1) * sizeof(double));

// 	/* Storing x(n) in a complex array to make it work with four1. 
// 	This is needed even though x(n) is purely real in this case. */
// 	for(i=0; i<Nx; i++)
// 	{
// 		X[2*i+1] = x[i];
// 		X[2*i+2] = 0.0;
// 	}
// 	 pad the remainder of the array with zeros (0 + 0 j) 
// 	for(i=Nx; i<NFFT; i++)
// 	{
// 		X[2*i+1] = 0.0;
// 		X[2*i+2] = 0.0;
// 	}

// 	printf("\nInput complex sequence (padded to next highest power of 2):\n");
// 	for(i=0; i<NFFT; i++)
// 	{
// 		printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2]);
// 	}

// 	/* calculate FFT */
// 	four1(X, NFFT, 1);

// 	printf("\nFFT:\n");
// 	for(i=0; i<NFFT; i++)
// 	{
// 		printf("X[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2]);
// 	}

// 	/* calculate IFFT */
// 	four1(X, NFFT, -1);

// 	/* normalize the IFFT */
// 	for(i=0; i<NFFT; i++)
// 	{
// 		X[2*i+1] /= NFFT;
// 		X[2*i+2] /= NFFT;
// 	}

// 	printf("\nComplex sequence reconstructed by IFFT:\n");
// 	for(i=0; i<NFFT; i++)
// 	{
// 		printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2]);
// 	}

// 	getchar();
// }

/*

Nx = 10
NFFT = 16

Input complex sequence (padded to next highest power of 2):
x[0] = (0.00 + j 0.00)
x[1] = (1.00 + j 0.00)
x[2] = (2.00 + j 0.00)
x[3] = (3.00 + j 0.00)
x[4] = (4.00 + j 0.00)
x[5] = (5.00 + j 0.00)
x[6] = (6.00 + j 0.00)
x[7] = (7.00 + j 0.00)
x[8] = (8.00 + j 0.00)
x[9] = (9.00 + j 0.00)
x[10] = (0.00 + j 0.00)
x[11] = (0.00 + j 0.00)
x[12] = (0.00 + j 0.00)
x[13] = (0.00 + j 0.00)
x[14] = (0.00 + j 0.00)
x[15] = (0.00 + j 0.00)

FFT:
X[0] = (45.00 + j 0.00)
X[1] = (-25.45 + j 16.67)
X[2] = (10.36 + j -3.29)
X[3] = (-9.06 + j -2.33)
X[4] = (4.00 + j 5.00)
X[5] = (-1.28 + j -5.64)
X[6] = (-2.36 + j 4.71)
X[7] = (3.80 + j -2.65)
X[8] = (-5.00 + j 0.00)
X[9] = (3.80 + j 2.65)
X[10] = (-2.36 + j -4.71)
X[11] = (-1.28 + j 5.64)
X[12] = (4.00 + j -5.00)
X[13] = (-9.06 + j 2.33)
X[14] = (10.36 + j 3.29)
X[15] = (-25.45 + j -16.67)

Complex sequence reconstructed by IFFT:
x[0] = (0.00 + j -0.00)
x[1] = (1.00 + j -0.00)
x[2] = (2.00 + j 0.00)
x[3] = (3.00 + j -0.00)
x[4] = (4.00 + j -0.00)
x[5] = (5.00 + j 0.00)
x[6] = (6.00 + j -0.00)
x[7] = (7.00 + j -0.00)
x[8] = (8.00 + j 0.00)
x[9] = (9.00 + j 0.00)
x[10] = (0.00 + j -0.00)
x[11] = (0.00 + j -0.00)
x[12] = (0.00 + j 0.00)
x[13] = (-0.00 + j -0.00)
x[14] = (0.00 + j 0.00)
x[15] = (0.00 + j 0.00)

*/









/////////////////////////////////////


void init_leds()
{
	int i;
	for (i = 0; i < 5; i++) {
		pinMode(rows[i], OUTPUT);
	}
	for (i = 0; i < 12; i++) {
		pinMode(cols[i], OUTPUT);
	}
}

void set_leds(int row, int col)
{
	int i;
	digitalWrite(cols[col], HIGH);
	for (i = 0; i < row; i++) {
		digitalWrite(rows[i], HIGH);
	}
}

void clear_leds()
{
	int i;
	for(i = 0; i < 5; i++)
	{
		digitalWrite(rows[i], LOW);
	}
	for(i = 0; i < 12; i++)
	{
		digitalWrite(cols[i], LOW);
	}
}

int rate = 44100;
int channels = 2;

int main(int argc, const char* argv[])
{

	wiringPiSetup();
	init_leds();

	int sock_fd, client_fd, err;
	struct sockaddr_in serv_addr, client_addr;
	int len = sizeof(client_addr);
	int received;

	serv_addr.sin_port = htons(7654);
	serv_addr.sin_family = AF_INET;

	int s;
	struct ifreq ifr = {};

	s = socket(PF_INET, SOCK_DGRAM, 0);

	strncpy(ifr.ifr_name, "eth0", sizeof(ifr.ifr_name));

	if (ioctl(s, SIOCGIFADDR, &ifr) >= 0)
		inet_pton(AF_INET, inet_ntoa(((struct sockaddr_in *)&ifr.ifr_addr)->sin_addr), &serv_addr.sin_addr.s_addr);
	else
		inet_pton(AF_INET, "192.168.1.9", &serv_addr.sin_addr.s_addr);

	/* SERVER DE MUZICA */

	int i;
	int rate = 44100;
	short buf[128];
	snd_pcm_t *playback_handle;
	snd_pcm_hw_params_t *hw_params;


	if ((err = snd_pcm_open (&playback_handle, "default", SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
		fprintf (stderr, "cannot open audio device %s (%s)\n", 
			 argv[1],
			 snd_strerror (err));
		exit (1);
	}
	   
	if ((err = snd_pcm_hw_params_malloc (&hw_params)) < 0) {
		fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n",
			 snd_strerror (err));
		exit (1);
	}
			 
	if ((err = snd_pcm_hw_params_any (playback_handle, hw_params)) < 0) {
		fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_set_access (playback_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
		fprintf (stderr, "cannot set access type (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_set_format (playback_handle, hw_params, SND_PCM_FORMAT_S16_LE)) < 0) {
		fprintf (stderr, "cannot set sample format (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_set_rate_near (playback_handle, hw_params, &rate, 0)) < 0) {
		fprintf (stderr, "cannot set sample rate (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	int channel = 2;

	if ((err = snd_pcm_hw_params_set_channels (playback_handle, hw_params, channel)) < 0) {
		fprintf (stderr, "cannot set channel count (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	snd_pcm_uframes_t frames;
	// frames = 32;
	int dir;
	// snd_pcm_hw_params_set_period_size_near(playback_handle, hw_params, &frames, &dir);

	if ((err = snd_pcm_hw_params (playback_handle, hw_params)) < 0) {
		fprintf (stderr, "cannot set parameters (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	snd_pcm_hw_params_get_period_size(hw_params, &frames, &dir);
	int buff_size = frames * 2 * channel; // size of buffer (1 period) in bytes
	char *buffer = (char *) malloc (buff_size);

	// printf("%dn", frames);
	
	snd_pcm_hw_params_free (hw_params);

	if ((err = snd_pcm_prepare (playback_handle)) < 0) {
		fprintf (stderr, "cannot prepare audio interface for use (%s)\n",
			 snd_strerror (err));
		exit (1);
	}

	/* END OF SERVER DE MUZICA */

	sock_fd = socket(AF_INET, SOCK_STREAM, 0);
	if(sock_fd < 0) {
		fprintf(stdout, "Cannot open socket.\n");
		return -1;
	}
	
	err = bind(sock_fd, (struct sockaddr*) &serv_addr, sizeof(serv_addr));
	if (err < 0) {
		fprintf(stdout, "Cannot bind socket to address.\n");
		return -1;
	}
	
	err = listen(sock_fd, 10);
	if (err < 0) {
		fprintf(stdout, "Cannot listen.\n");
		return -1;
	}

	printf("Listening...\n");

	client_fd = accept(sock_fd, (struct sockaddr*) &client_addr, &len);

	if (client_fd < 0) {
		fprintf(stdout, "Could not connect to the client!\n");
		return -1;
	}

	printf("New client connected!\n");
	double *samples = (double *) malloc ((2 * frames + 1) * sizeof(double));
	float levels[12];

	while (1) {
		received = read(client_fd, buffer, buff_size);

		if (received > 0) {
			// if (received == frames / 4)
			//printf("Received: %d\n", received);
			if (frames == -EPIPE) {
				snd_pcm_prepare(pcm.handle);
				if ((err = snd_pcm_writei(playback_handle, buffer, frames)) <= 0) {
					fprintf(stderr, "write to audio interface failed (%s)\n", snd_strerror(err));
				}
			}

			
			// printf("Play\n");

			/* create input vector for fft */
			/* samples are real, imaginary part is zero */
			for (i = 0; i < received / 4; i++) {
				samples[2 * i + 1] = buffer[2 * i] / 32768.0f; // divide by ... to convert from short to float
				samples[2 * i + 2] = 0.0;
			}

			four1(samples, received / 4, 1);
			float max = 0;

			for (i = 0; i < 12; i++) {
				int j;
				int start = i * (received / 4 / 12);
				for (j = start; j < start + (received / 4/ 12); j++) {
					levels[i] = sqrt(samples[2 * j + 1] * samples[2 * j + 1] + samples[2 * j + 2] * samples[2 * j + 2]);
				}
				// levels[i] /= frames / 12;
				if (levels[i] > max)
					max = levels[i];
			}

			for (i = 0; i < 12; i++) {
				levels[i] /= max;
				levels[i] *= 5;
				
				int lvl = (int) levels[i]; /*NOT FUCKING WORKING */
				clear_leds();
				set_leds(lvl, i); // row, col
				//
				printf("Column: %d, Level: %d\n", i, lvl);
			}
		}
		else
		{
			clear_leds();
		}
	}

	snd_pcm_close (playback_handle);
	return 0;
}




// to copy to raspberry: scp [-r] path/to/file user@host:path/to/remote/folder
// to copy from raspberry: scp [-r] user@host:/path/to/remote/file /path/to/folder
