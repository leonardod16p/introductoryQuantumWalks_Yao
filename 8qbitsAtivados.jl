# ======================================================================
# CAMINHADA QUÂNTICA DTW - 8 QUBITS
# Usando o "motor" XY validado do teste de 2 qubits
# ======================================================================

using Yao
using CairoMakie
using LinearAlgebra

println("=== CAMINHADA QUÂNTICA DTW - 8 QUBITS ===")
println("Motor XY validado: J(X₁X₂ + Y₁Y₂) para troca real")
println("Sistema: 8 qubits em cadeia linear")
println("Protocolo: Alternância W₀ (intra-dímero) ↔ W₁ (inter-dímero)")

# --- 1. CONFIGURAÇÃO DO SISTEMA ---
N = 8  # 8 qubits
println("✓ Sistema de $N qubits criado")

# Estado inicial: excitação no primeiro qubit
initial_state = product_state(bit"00000001")
println("✓ Estado inicial: |10000000⟩ (primeiro qubit excitado)")

# --- 2. FUNÇÃO DO MOTOR XY VALIDADO ---
function create_xy_coupling(N, i, j, J_coupling)
    """
    Motor XY entre qubits i e j: H = J * (XᵢXⱼ + YᵢYⱼ)
    Esta é a física validada que permite troca real |10⟩ ↔ |01⟩
    """
    Xi = put(N, i => X)
    Xj = put(N, j => X)
    Yi = put(N, i => Y)
    Yj = put(N, j => Y)
    
    return J_coupling * (Xi * Xj + Yi * Yj)
end

# --- 3. HAMILTONIANOS DOS PASSOS DTW ---

function create_step_hamiltonian(N, step_type, J_coupling)
    """
    Cria Hamiltoniano para um passo da DTW
    step_type: "W0" (intra-dímero) ou "W1" (inter-dímero)
    """
    
    if step_type == "W0"
        # W₀: Acoplamentos intra-dímero (1-2), (3-4), (5-6), (7-8)
        pairs = [(1,2), (3,4), (5,6), (7,8)]
        println("  Construindo W₀: acoplamentos intra-dímero")
    elseif step_type == "W1"
        # W₁: Acoplamentos inter-dímero (2-3), (4-5), (6-7)
        pairs = [(2,3), (4,5), (6,7)]
        println("  Construindo W₁: acoplamentos inter-dímero")
    else
        error("step_type deve ser 'W0' ou 'W1'")
    end
    
    # ======================= A CORREÇÃO FINAL ESTÁ AQUI ==========================
    # Usamos a função `sum` do Julia para somar os Hamiltonianos de todos os pares
    # de uma só vez. É mais limpo, mais eficiente e evita o erro de repetição.
    H_total = sum(create_xy_coupling(N, i, j, J_coupling) for (i, j) in pairs)
    
    return H_total
end

# Parâmetros físicos (mesmos do teste validado)
J = 2π * 1.0e6  # Hz - força do acoplamento
τ = π / (4*J)  # Tempo de cada passo (π/2 pulse)

println("\n✓ Parâmetros físicos:")
println("  Acoplamento J: $(J/1e6) MHz")
println("  Tempo por passo: $(τ*1e6) μs")

# Criar os Hamiltonianos
H_W0 = create_step_hamiltonian(N, "W0", J)
H_W1 = create_step_hamiltonian(N, "W1", J)

println("✓ Hamiltonianos W₀ e W₁ construídos")

# --- 4. FUNÇÃO PARA CALCULAR PROBABILIDADES ---
function calculate_all_probabilities(state, N)
    """Calcula probabilidade de cada qubit estar em |1⟩"""
    probs = zeros(Float64, N)
    
    for i in 1:N
        z_op = put(N, i => Z)
        z_exp = real(expect(z_op, state))
        probs[i] = (1 - z_exp) / 2  # Converter ⟨Z⟩ → P(|1⟩)
    end
    
    return probs
end

# --- 5. SIMULAÇÃO DA CAMINHADA QUÂNTICA ---

println("\n🔬 Iniciando simulação da caminhada quântica DTW...")

# Parâmetros da simulação
n_steps = 16  # Número de passos DTW
n_points_per_step = 50  # Resolução temporal dentro de cada passo

# Arrays para armazenar resultados
time_points = Float64[]
prob_evolution = Vector{Float64}[]  # Array de arrays para cada qubit
step_markers = Float64[]  # Marcar início de cada passo
step_labels = String[]    # Labels dos passos

# Estado atual
current_state = copy(initial_state)

# Tempo atual
current_time = 0.0

println("  Configuração: $n_steps passos, $n_points_per_step pontos por passo")

for step in 1:n_steps
    # Determinar tipo de passo
    step_type = (step % 2 == 1) ? "W0" : "W1"
    H_step = (step_type == "W0") ? H_W0 : H_W1
    
    println("  Passo $step: $step_type")
    push!(step_markers, current_time)
    push!(step_labels, step_type)
    
    # Evolução temporal dentro do passo
    dt = τ / (n_points_per_step - 1)
    
    for i in 1:n_points_per_step
        # Calcular probabilidades atuais
        probs = calculate_all_probabilities(current_state, N)
        
        # Armazenar resultados
        push!(time_points, current_time)
        push!(prob_evolution, copy(probs))
        
        # Evoluir para próximo ponto (exceto no último)
        if i < n_points_per_step
            U_step = time_evolve(H_step, dt)
# ======================= A CORREÇÃO ESTÁ AQUI ==========================
            # Usamos `global` para dizer ao Julia para modificar as variáveis
            # que foram criadas fora do loop.
            global current_state = apply!(current_state, U_step)
            global current_time += dt
            # ==============================================
        end
    end
    
    # Debug a cada 4 passos
    if step % 4 == 0
        final_probs = prob_evolution[end]
        max_prob = maximum(final_probs)
        max_qubit = argmax(final_probs)
        println("    → t=$(round(current_time*1e6, digits=1))μs: Máx P$max_qubit = $(round(max_prob, digits=3))")
    end
end

println("✓ Simulação concluída!")

# --- 6. ANÁLISE DOS RESULTADOS ---

# Converter para matriz para facilitar análise
prob_matrix = hcat(prob_evolution...)  # Cada coluna é um tempo, cada linha um qubit
time_μs = time_points .* 1e6

println("\n" * "="^60)
println("📊 ANÁLISE DA CAMINHADA QUÂNTICA:")
println("="^60)

# Calcular métricas de propagação
for qubit in 1:N
    max_prob = maximum(prob_matrix[qubit, :])
    max_time_idx = argmax(prob_matrix[qubit, :])
    max_time = time_μs[max_time_idx]
    
    println("Qubit $qubit: Máx probabilidade = $(round(max_prob, digits=3)) em t=$(round(max_time, digits=1))μs")
end

# Verificar conservação
total_probs = sum(prob_matrix, dims=1)[:]  # Soma por tempo
avg_conservation = sum(total_probs) / length(total_probs)
println("\nConservação média: $(round(avg_conservation, digits=6))")

# Métricas de propagação
initial_center = sum([i * prob_matrix[i, 1] for i in 1:N])
final_center = sum([i * prob_matrix[i, end] for i in 1:N])
propagation_distance = final_center - initial_center

println("Posição média inicial: $(round(initial_center, digits=2))")
println("Posição média final: $(round(final_center, digits=2))")
println("Distância de propagação: $(round(propagation_distance, digits=2)) qubits")

if propagation_distance > 2.0
    println("✅ CAMINHADA DETECTADA! A excitação se propagou pela cadeia")
else
    println("⚠️  Propagação limitada. A excitação ainda está localizada")
end

println("="^60)

# --- 7. VISUALIZAÇÃO COMPLETA ---

fig = Figure(size = (1400, 1000))

# Gráfico principal: Evolução de todas as probabilidades
ax1 = Axis(fig[1, 1:2],
    title = "Quantum Walk DTW - 8 QUBITS\n |10000000⟩ ↔ |00000001⟩",
    xlabel = "Tempo (μs)",
    ylabel = "Probabilidade P(|1⟩)",
    limits = (nothing, nothing, -0.05, 1.05)
)

# Cores para cada qubit
colors = [:red, :blue, :green, :orange, :purple, :brown, :pink, :gray]

# Plotar evolução de cada qubit
for qubit in 1:N
    lines!(ax1, time_μs, prob_matrix[qubit, :], 
           linewidth=3, color=colors[qubit], label="Qubit $qubit")
end

# Marcar transições entre passos
for (i, t_marker) in enumerate(step_markers[2:end])  # Pular primeiro ponto
    vlines!(ax1, [t_marker], color=:gray, linestyle=:dash, alpha=0.6)
    if i % 2 == 1
        text!(ax1, t_marker, 0.9, text="W₀", align=(:center, :bottom), fontsize=10)
    else
        text!(ax1, t_marker, 0.85, text="W₁", align=(:center, :bottom), fontsize=10)
    end
end

axislegend(ax1, position=:rt, nbanks=2)

# Gráfico de conservação
ax2 = Axis(fig[2, 1],
    title = "Conservação Total",
    xlabel = "Tempo (μs)",
    ylabel = "∑ Pᵢ",
    limits = (nothing, nothing, 0.98, 1.02)
)

lines!(ax2, time_μs, total_probs, linewidth=3, color=:purple)
hlines!(ax2, [1.0], color=:black, linestyle=:solid, alpha=0.7)

# Gráfico da posição média (centro de massa)
ax3 = Axis(fig[2, 2],
    title = "Posição Média da Excitação",
    xlabel = "Tempo (μs)",
    ylabel = "⟨Posição⟩",
    limits = (nothing, nothing, 0.5, 8.5)
)

# Calcular posição média ao longo do tempo
mean_positions = [sum([i * prob_matrix[i, t] for i in 1:N]) for t in 1:length(time_points)]
lines!(ax3, time_μs, mean_positions, linewidth=3, color=:darkgreen)

# Marcar posições dos qubits
hlines!(ax3, collect(1:N), color=:gray, linestyle=:dot, alpha=0.3)

println("✓ Gráficos criados")

# --- 8. DIAGNÓSTICO FINAL ---

println("\n" * "="^60)
println("🏁 DIAGNÓSTICO FINAL DA CAMINHADA QUÂNTICA:")
println("="^60)

# Verificar se houve propagação efetiva
qubits_activated = sum([maximum(prob_matrix[i, :]) > 0.1 for i in 1:N])
propagation_efficiency = (final_center - 1) / (N - 1)  # Normalizado [0,1]

println("🎯 Métricas de Propagação:")
println("  Qubits significativamente ativados: $qubits_activated de $N")
println("  Eficiência de propagação: $(round(propagation_efficiency*100, digits=1))%")
println("  Distância percorrida: $(round(propagation_distance, digits=1)) posições")

println("\n🔬 Interpretação Física:")
if qubits_activated >= 4 && propagation_efficiency > 0.3
    println("✅ CAMINHADA QUÂNTICA REALIZADA COM SUCESSO!")
    println("   → A excitação se propagou através da cadeia")
    println("   → O protocolo DTW funciona com o motor XY")
    println("   → Múltiplos qubits foram ativados sequencialmente")
elseif qubits_activated >= 3
    println("⚠️  PROPAGAÇÃO PARCIAL DETECTADA")
    println("   → Alguns qubits foram ativados mas propagação limitada")
    println("   → Pode precisar ajustar parâmetros temporais")
else
    println("❌ PROPAGAÇÃO INSUFICIENTE")
    println("   → A excitação permanece muito localizada")
    println("   → Verificar acoplamentos ou timing")
end

println("\n🎲 Validação da Física:")
if abs(avg_conservation - 1.0) < 0.01
    println("✅ Conservação de probabilidade mantida")
else
    println("⚠️  Erro de conservação: $(abs(avg_conservation - 1.0))")
end

println("\n" * "="^60)
println("🎉 CONCLUSÃO:")

if qubits_activated >= 4
    println("🏆 MISSÃO CUMPRIDA!")
    println("   • Motor XY validado funcionou perfeitamente")
    println("   • Caminhada quântica DTW implementada com sucesso") 
    println("   • A excitação 'caminhou' pela cadeia de qubits")
    println("   • Física quântica conservativa mantida")
else
    println("🔧 AJUSTES NECESSÁRIOS:")
    println("   • Motor XY funciona mas protocolo DTW precisa otimização")
    println("   • Considerar ajustar tempos τ ou número de passos")
    println("   • A física está correta mas parâmetros podem melhorar")
end
println("="^60)

display(fig)

println("\n📊 OBSERVAÇÕES NO GRÁFICO:")
println("   • Picos se movem da esquerda para direita? = Caminhada funcionando")
println("   • Excitação fica presa no início? = Parâmetros a ajustar")
println("   • Conservação = 1? = Física correta")
println("   • Posição média cresce? = Propagação efetiva")

println("\n🧪 Este é o teste definitivo da caminhada quântica DTW!")
println("💡 Motor XY do teste de 2 qubits agora expandido para 8 qubits!")