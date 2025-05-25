// グローバルスコープにチャートオブジェクトを保持
let freqChart = null;
let zPlaneChart = null;
let impulseChart = null;
// Plotlyの3Dグラフは newPlot で都度生成するので、オブジェクト保持は必須ではない

const initialCoefficients = {
    b0: 1.0, b1: 0.0, b2: 0.0,
    a0: 1.0, a1: 0.0, a2: 0.0
};

document.addEventListener('DOMContentLoaded', () => {
    // 入力フィールドにイベントリスナーを設定
    const coeffIds = ['b0', 'b1', 'b2', 'a1', 'a2']; // a0 は固定
    coeffIds.forEach(id => {
        const element = document.getElementById(id);
        if (element) {
            // 'input' イベントで頻繁に更新、デバウンス処理を追加
            let debounceTimer;
            element.addEventListener('input', () => {
                clearTimeout(debounceTimer);
                debounceTimer = setTimeout(updateAllPlots, 250); // 250msのデバウンス
            });
        }
    });

    // リセットボタンの処理
    const resetButton = document.getElementById('resetCoeffs');
    if (resetButton) {
        resetButton.addEventListener('click', () => {
            document.getElementById('b0').value = initialCoefficients.b0;
            document.getElementById('b1').value = initialCoefficients.b1;
            document.getElementById('b2').value = initialCoefficients.b2;
            document.getElementById('a1').value = initialCoefficients.a1;
            document.getElementById('a2').value = initialCoefficients.a2;
            updateAllPlots();
        });
    }

    // 初期描画
    setInitialCoefficients();
    updateAllPlots();
});

function setInitialCoefficients() {
    document.getElementById('b0').value = initialCoefficients.b0;
    document.getElementById('b1').value = initialCoefficients.b1;
    document.getElementById('b2').value = initialCoefficients.b2;
    // a0 is disabled and fixed to 1.0
    document.getElementById('a1').value = initialCoefficients.a1;
    document.getElementById('a2').value = initialCoefficients.a2;
}


// フィルタ係数を取得
function getCoefficients() {
    const b0 = parseFloat(document.getElementById('b0').value) || 0;
    const b1 = parseFloat(document.getElementById('b1').value) || 0;
    const b2 = parseFloat(document.getElementById('b2').value) || 0;
    const a0 = 1.0; // a0 は 1 に正規化されていると仮定
    const a1 = parseFloat(document.getElementById('a1').value) || 0;
    const a2 = parseFloat(document.getElementById('a2').value) || 0;
    return { b0, b1, b2, a0, a1, a2 };
}

// --- 計算関数 ---

// 周波数特性の計算
function calculateFrequencyResponse(coeffs, numPoints = 256) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    const magnitudes = new Array(numPoints);
    const phases = new Array(numPoints);
    const frequencies = new Array(numPoints);

    for (let i = 0; i < numPoints; i++) {
        const omega = (Math.PI * i) / (numPoints - 1); // 0 から π
        frequencies[i] = omega / Math.PI; // 正規化周波数 (0 から 1, fs/2に対応)

        const cosOmega = Math.cos(omega);
        const sinOmega = Math.sin(omega);
        const cos2Omega = Math.cos(2 * omega);
        const sin2Omega = Math.sin(2 * omega);

        // 分子: N(z) = b0 + b1*z^-1 + b2*z^-2
        const numReal = b0 + b1 * cosOmega + b2 * cos2Omega;
        const numImag = -b1 * sinOmega - b2 * sin2Omega;

        // 分母: D(z) = a0 + a1*z^-1 + a2*z^-2
        const denReal = a0 + a1 * cosOmega + a2 * cos2Omega;
        const denImag = -a1 * sinOmega - a2 * sin2Omega;

        const denMagnitudeSq = denReal * denReal + denImag * denImag;

        if (denMagnitudeSq < 1e-12) { // 極に近い場合、非常に大きな値
            magnitudes[i] = 1e6; // 大きな値 (dBにすると120dB)
            phases[i] = (numImag * denReal - numReal * denImag === 0) ? 0 : Math.atan2(numImag * denReal - numReal * denImag, numReal * denReal + numImag * denImag);
        } else {
            const hReal = (numReal * denReal + numImag * denImag) / denMagnitudeSq;
            const hImag = (numImag * denReal - numReal * denImag) / denMagnitudeSq;
            magnitudes[i] = Math.sqrt(hReal * hReal + hImag * hImag);
            phases[i] = Math.atan2(hImag, hReal);
        }
    }
    return { frequencies, magnitudes, phases };
}

// 二次方程式の解 (az^2 + bz + c = 0)
function solveQuadratic(a, b, c) {
    if (Math.abs(a) < 1e-9) { // ほぼ0の場合、一次方程式
        if (Math.abs(b) < 1e-9) return []; // 解なしまたは不定
        return [{ real: -c / b, imag: 0 }];
    }
    const discriminant = b * b - 4 * a * c;
    const twoA = 2 * a;
    const roots = [];
    if (discriminant >= 0) {
        roots.push({ real: (-b + Math.sqrt(discriminant)) / twoA, imag: 0 });
        roots.push({ real: (-b - Math.sqrt(discriminant)) / twoA, imag: 0 });
    } else {
        roots.push({ real: -b / twoA, imag: Math.sqrt(-discriminant) / twoA });
        roots.push({ real: -b / twoA, imag: -Math.sqrt(-discriminant) / twoA });
    }
    return roots;
}

// 極と零点の計算
function calculatePolesAndZeros(coeffs) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    // 零点: b0*z^2 + b1*z + b2 = 0 (z^-1の多項式 b0 + b1*x + b2*x^2 = 0, x=z^-1 の根を求め、z=1/x)
    // または、H(z) の分子を z^2 で割って b0*z^2 + b1*z + b2 = 0 の根を求める
    // 伝達関数 H(z) = (b0 + b1 z^-1 + b2 z^-2) / (a0 + a1 z^-1 + a2 z^-2)
    // 分子 = 0: b0 z^2 + b1 z + b2 = 0 (両辺に z^2 をかける)
    // 分母 = 0: a0 z^2 + a1 z + a2 = 0 (両辺に z^2 をかける)
    const zeros = solveQuadratic(b0, b1, b2);
    const poles = solveQuadratic(a0, a1, a2); // a0=1
    return { zeros, poles };
}


// インパルス応答の計算
function calculateImpulseResponse(coeffs, length = 64) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs; // a0 は 1 と仮定
    const y = new Array(length).fill(0);
    const x = new Array(length).fill(0);
    x[0] = 1; // インパルス入力

    for (let n = 0; n < length; n++) {
        let val = b0 * x[n];
        if (n >= 1) val += b1 * x[n - 1];
        if (n >= 2) val += b2 * x[n - 2];

        if (n >= 1) val -= a1 * y[n - 1];
        if (n >= 2) val -= a2 * y[n - 2];
        
        y[n] = val / a0; // a0で正規化 (a0=1なので実質そのまま)
    }
    return y;
}

// Z平面3Dデータ計算
function calculateZPlane3DData(coeffs, gridSize = 25) { // gridSizeを少し小さめに
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    const xGridValues = []; // 実数部のグリッド値
    const yGridValues = []; // 虚数部のグリッド値
    const zMagnitudes = []; // 各グリッド点での振幅(dB)

    const range = 1.5; // z平面の表示範囲 (-range から +range)

    for (let i = 0; i < gridSize; i++) {
        const re = -range + (2 * range * i) / (gridSize - 1);
        xGridValues.push(re);
        yGridValues.push(re); // yGridValuesは実際には虚数部のimに対応するが、Plotlyではx,yとして渡す
    }

    for (let j = 0; j < gridSize; j++) { // y-axis (imaginary part)
        const im = yGridValues[j]; // -range + (2 * range * j) / (gridSize - 1);
        const rowMagnitudes = [];
        for (let i = 0; i < gridSize; i++) { // x-axis (real part)
            const re = xGridValues[i];

            const modSqZ = re * re + im * im;

            let numReal, numImag, denReal, denImag;

            if (modSqZ < 1e-9) { // zがほぼ0の場合
                // z^-1, z^-2 が発散する可能性。H(z)の定義から考える。
                // H(z) = (b0 z^2 + b1 z + b2) / (a0 z^2 + a1 z + a2) * z^0 (if a0, b0 are for z^0 term)
                // H(z) = (b0 + b1 z^-1 + b2 z^-2) / (a0 + a1 z^-1 + a2 z^-2)
                // z->0 のとき、z^-1, z^-2 の項が支配的。
                // もし b2, a2 が非ゼロなら H(z) -> b2/a2
                // もし b2=0, a2=0 で b1, a1 が非ゼロなら H(z) -> (b1 z^-1)/(a1 z^-1) = b1/a1
                // もし b2=b1=0, a2=a1=0 なら H(z) -> b0/a0
                if (Math.abs(a2) > 1e-9) {
                    numReal = b2; numImag = 0; denReal = a2; denImag = 0;
                } else if (Math.abs(a1) > 1e-9) {
                    numReal = b1; numImag = 0; denReal = a1; denImag = 0;
                } else if (Math.abs(a0) > 1e-9) {
                    numReal = b0; numImag = 0; denReal = a0; denImag = 0;
                } else { // 分母が0になるケース
                     rowMagnitudes.push(60); // 大きな値 (dB)
                    continue;
                }

            } else {
                const zInvRe = re / modSqZ;
                const zInvIm = -im / modSqZ;

                const zInvSqRe = (re * re - im * im) / (modSqZ * modSqZ);
                const zInvSqIm = (-2 * re * im) / (modSqZ * modSqZ);
                
                numReal = b0 + b1 * zInvRe + b2 * zInvSqRe;
                numImag = b1 * zInvIm + b2 * zInvSqIm;

                denReal = a0 + a1 * zInvRe + a2 * zInvSqRe;
                denImag = a1 * zInvIm + a2 * zInvSqIm;
            }
            
            const denMagSq = denReal * denReal + denImag * denImag;

            if (denMagSq < 1e-12) {
                rowMagnitudes.push(60); // 極に近い場合、大きな値 (dB)
            } else {
                const hReal = (numReal * denReal + numImag * denImag) / denMagSq;
                const hImag = (numImag * denReal - numReal * denImag) / denMagSq;
                const magnitude = Math.sqrt(hReal * hReal + hImag * hImag);
                rowMagnitudes.push(Math.min(60, Math.max(-60, 20 * Math.log10(Math.max(magnitude, 1e-6))))); // dBスケール, 上下限設定
            }
        }
        zMagnitudes.push(rowMagnitudes);
    }
    return { x: xGridValues, y: yGridValues, z: zMagnitudes };
}


// --- 描画関数 ---

// 周波数特性グラフ
function plotFrequencyResponse(data) {
    const { frequencies, magnitudes, phases } = data;
    const ctx = document.getElementById('freqResponseChart').getContext('2d');

    const magDb = magnitudes.map(m => Math.min(60, Math.max(-100, 20 * Math.log10(Math.max(m, 1e-5))))); // dB, 上下限
    const phaseDeg = phases.map(p => p * 180 / Math.PI);

    if (freqChart) {
        freqChart.data.labels = frequencies;
        freqChart.data.datasets[0].data = magDb;
        freqChart.data.datasets[1].data = phaseDeg;
        freqChart.update();
    } else {
        freqChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: frequencies,
                datasets: [
                    {
                        label: '振幅特性 (dB)',
                        data: magDb,
                        borderColor: 'rgb(75, 192, 192)',
                        yAxisID: 'y-axis-mag',
                        tension: 0.1,
                        pointRadius: 0,
                    },
                    {
                        label: '位相特性 (度)',
                        data: phaseDeg,
                        borderColor: 'rgb(255, 99, 132)',
                        yAxisID: 'y-axis-phase',
                        tension: 0.1,
                        pointRadius: 0,
                    }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: { title: { display: true, text: '正規化周波数 (×π rad/sample)' } },
                    'y-axis-mag': { type: 'linear', position: 'left', title: { display: true, text: '振幅 (dB)' } },
                    'y-axis-phase': { type: 'linear', position: 'right', title: { display: true, text: '位相 (度)' }, grid: { drawOnChartArea: false } }
                },
                interaction: { mode: 'index', intersect: false },
            }
        });
    }
}

// z平面 極・零点グラフ
function plotZPlanePolesZeros(data) {
    const { poles, zeros } = data;
    const ctx = document.getElementById('zPlaneChart').getContext('2d');

    const poleData = poles.map(p => ({ x: p.real, y: p.imag }));
    const zeroData = zeros.map(z => ({ x: z.real, y: z.imag }));

    // 単位円のデータ
    const unitCircle = [];
    for (let i = 0; i <= 360; i += 5) {
        const angle = i * Math.PI / 180;
        unitCircle.push({ x: Math.cos(angle), y: Math.sin(angle) });
    }

    if (zPlaneChart) {
        zPlaneChart.data.datasets[0].data = unitCircle;
        zPlaneChart.data.datasets[1].data = poleData;
        zPlaneChart.data.datasets[2].data = zeroData;
        zPlaneChart.update();
    } else {
        zPlaneChart = new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [
                    {
                        label: '単位円',
                        data: unitCircle,
                        borderColor: 'rgba(150, 150, 150, 0.8)',
                        borderWidth: 1.5,
                        pointRadius: 0,
                        showLine: true,
                        tension: 0.1,
                        fill: false,
                        order: 3 // 描画順序（奥）
                    },
                    {
                        label: '極 (Poles)',
                        data: poleData,
                        backgroundColor: 'rgba(255, 0, 0, 0.7)',
                        pointStyle: 'crossRot', // 'X' marks
                        radius: 8,
                        hoverRadius: 10,
                        order: 1
                    },
                    {
                        label: '零点 (Zeros)',
                        data: zeroData,
                        backgroundColor: 'rgba(0, 0, 255, 0.7)',
                        pointStyle: 'circle', // 'O' marks
                        radius: 6,
                        hoverRadius: 8,
                        borderColor: 'rgba(0,0,255,1)',
                        borderWidth: 1,
                        order: 2
                    }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false, // アスペクト比をコンテナに合わせる
                aspectRatio: 1, // アスペクト比を1:1に
                scales: {
                    x: {
                        title: { display: true, text: '実数部 (Re)' },
                        min: -1.5, max: 1.5,
                        grid: { color: '#ddd' }
                    },
                    y: {
                        title: { display: true, text: '虚数部 (Im)' },
                        min: -1.5, max: 1.5,
                        grid: { color: '#ddd' }
                    }
                },
                plugins: {
                    legend: { display: true }
                }
            }
        });
    }
}

// インパルス応答グラフ
function plotImpulseResponse(impulseData) {
    const ctx = document.getElementById('impulseResponseChart').getContext('2d');
    const labels = impulseData.map((_, i) => i);

    if (impulseChart) {
        impulseChart.data.labels = labels;
        impulseChart.data.datasets[0].data = impulseData;
        impulseChart.update();
    } else {
        impulseChart = new Chart(ctx, {
            type: 'bar', // または 'line'
            data: {
                labels: labels,
                datasets: [{
                    label: 'インパルス応答 y[n]',
                    data: impulseData,
                    backgroundColor: 'rgba(54, 162, 235, 0.6)',
                    borderColor: 'rgb(54, 162, 235)',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: { title: { display: true, text: 'サンプル (n)' } },
                    y: { title: { display: true, text: '振幅' } }
                }
            }
        });
    }
}

// Z平面3D振幅特性グラフ
async function plotZPlane3D(data) {
    const { x, y, z } = data;
    const plotlyLoader = document.getElementById('plotlyLoader');
    if (plotlyLoader) plotlyLoader.style.display = 'block'; // ローダー表示

    // Plotlyの処理は非同期に感じる場合があるため、少し遅延させる
    await new Promise(resolve => setTimeout(resolve, 50));


    const plotData = [{
        type: 'surface',
        x: x,
        y: y,
        z: z,
        colorscale: 'Viridis',
        contours: {
            z: { show: true, usecolormap: true, project: { z: true } }
        },
        showscale: true, // カラーバーを表示
        colorbar: {title: '振幅 (dB)'}
    }];

    const layout = {
        title: 'Z平面 3D振幅特性',
        scene: {
            xaxis: { title: '実数部 (Re)', range: [-1.5, 1.5] },
            yaxis: { title: '虚数部 (Im)', range: [-1.5, 1.5] },
            zaxis: { title: '振幅 (dB)', range: [-60, 60] }, // Z軸の範囲を固定
            camera: { eye: { x: 1.5, y: 1.5, z: 1.2 } },
            aspectmode: 'cube' // アスペクト比を立方体に
        },
        autosize: true,
        margin: { l: 40, r: 20, b: 40, t: 60 }
    };
    
    try {
        await Plotly.newPlot('zPlane3DChart', plotData, layout, {responsive: true});
    } catch (error) {
        console.error("Plotly 3D plot error:", error);
        const chartDiv = document.getElementById('zPlane3DChart');
        if(chartDiv) chartDiv.innerHTML = "<p style='color:red; text-align:center;'>3Dグラフの描画に失敗しました。</p>";
    } finally {
        if (plotlyLoader) plotlyLoader.style.display = 'none'; // ローダー非表示
    }
}


// --- 全プロット更新 ---
let isPlotting3D = false; // 3Dプロット中のフラグ

async function updateAllPlots() {
    const coeffs = getCoefficients();

    // 周波数特性
    const freqData = calculateFrequencyResponse(coeffs);
    plotFrequencyResponse(freqData);

    // 極・零点
    const pzData = calculatePolesAndZeros(coeffs);
    plotZPlanePolesZeros(pzData);

    // インパルス応答
    const impulseData = calculateImpulseResponse(coeffs);
    plotImpulseResponse(impulseData);

    // Z平面3D (非同期処理で、多重実行を避ける)
    if (!isPlotting3D) {
        isPlotting3D = true;
        try {
            const zPlane3DData = calculateZPlane3DData(coeffs);
            await plotZPlane3D(zPlane3DData);
        } catch (error) {
            console.error("Error during 3D plot update:", error);
        } finally {
            isPlotting3D = false;
        }
    }
}

