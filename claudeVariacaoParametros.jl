# ======================================================================
# CAMINHADA QUÂNTICA DTW - PULSOS VARIÁVEIS NO TEMPO
# Parâmetros J(t) e τ(t) que evoluem durante a simulação
# ======================================================================

using Yao
using CairoMakie
using LinearAlgebra

println("=== CAMINHADA QUÂNTICA DTW - PULSOS VARIÁVEIS NO TEMPO ===")
println("Motor XY: J(t)(X₁X₂ + Y₁Y₂) com parâmetros dinâmicos")
println("Sistema: 8 qubits em cadeia linear")
println("NOVIDADE: J(t) e τ(t) variam durante a evolução!")

# --- 1. CONFIGURAÇÃO DO SISTEMA ---
N = 8
println("✓ Sistema de $N qubits criado")

initial_state = product_state(bit"00000001")
println("✓ Estado inicial: |10000000⟩ (primeiro qubit excitado)")

# --- 2. FUNÇÕES DE VARIAÇÃO TEMPORAL ---

# Perfis de variação temporal disponíveis
function constant_profile(t, base_value, kwargs...)
    """Perfil constante (referência)"""
    return base_value
end

function linear_ramp(t, base_value; t_max=1e-6, final_factor=2.0)
    """Rampa linear: cresce linearmente até final_factor"""
    factor = 1.0 + (final_factor - 1.0) * min(t / t_max, 1.0)
    return base_value * factor
end

function exponential_decay(t, base_value; decay_rate=1e6)
    """Decaimento exponencial"""
    return base_value * exp(-decay_rate * t)
end

function sinusoidal_modulation(t, base_value; freq=5e6, amplitude=0.5)
    """Modulação senoidal"""
    return base_value * (1.0 + amplitude * sin(2π * freq * t))
end

function step_function(t, base_value; t_switch=5e-7, factor_before=1.0, factor_after=3.0)
    """Função degrau: muda abruptamente em t_switch"""
    factor = t < t_switch ? factor_before : factor_after
    return base_value * factor
end

function gaussian_pulse(t, base_value; t_center=5e-7, width=1e-7, amplitude=2.0)
    """Pulso gaussiano centrado em t_center"""
    gaussian = amplitude * exp(-((t - t_center) / width)^2)
    return base_value * (1.0 + gaussian)
end

function chirp_modulation(t, base_value; freq_start=1e6, freq_end=10e6, t_max=1e-6)
    """Chirp: frequência varia linearmente (útil para controle adiabático)"""
    freq_t = freq_start + (freq_end - freq_start) * min(t / t_max, 1.0)
    return base_value * (1.0 + 0.3 * sin(2π * freq_t * t))
end

# --- 3. MOTOR XY E HAMILTONIANOS ---
function create_xy_coupling(N, i, j, J_coupling)
    """Motor XY entre qubits i e j"""
    Xi = put(N, i => X)
    Xj = put(N, j => X)
    Yi = put(N, i => Y)  
    Yj = put(N, j => Y)
    return J_coupling * (Xi * Xj + Yi * Yj)
end

function create_step_hamiltonian(N, step_type, J_coupling)
    """Cria Hamiltoniano para um passo da DTW"""
    if step_type == "W0"
        pairs = [(1,2), (3,4), (5,6), (7,8)]  # intra-dímero
    elseif step_type == "W1"
        pairs = [(2,3), (4,5), (6,7)]  # inter-dímero
    else
        error("step_type deve ser 'W0' ou 'W1'")
    end
    
    return sum(create_xy_coupling(N, i, j, J_coupling) for (i, j) in pairs)
end

# --- 4. FUNÇÃO DE SIMULAÇÃO COM PULSOS VARIÁVEIS ---

function simulate_variable_quantum_walk(
    J_profile_func, τ_profile_func, 
    J_base=1.0e6, τ_base_factor=10.0,
    n_steps=32, n_points_per_step=25;
    J_kwargs=NamedTuple(), τ_kwargs=NamedTuple()
)
    """
    Simula caminhada quântica com parâmetros variáveis no tempo
    
    Args:
        J_profile_func: função para J(t) 
        τ_profile_func: função para τ(t)
        J_base: valor base para J (Hz)
        τ_base_factor: fator base para τ
        J_kwargs, τ_kwargs: parâmetros específicos para cada perfil
    """
    
    println("\n" * "="^50)
    println("🚀 SIMULANDO PULSOS VARIÁVEIS NO TEMPO")
    println("="^50)
    println("Perfil J(t): $(J_profile_func)")
    println("Perfil τ(t): $(τ_profile_func)")
    println("Parâmetros J: $J_kwargs")
    println("Parâmetros τ: $τ_kwargs")
    
    # Arrays para resultados
    time_points = Float64[]
    prob_evolution = Vector{Float64}[]
    J_evolution = Float64[]
    τ_evolution = Float64[]
    step_markers = Float64[]
    step_labels = String[]
    
    current_state = copy(initial_state)
    current_time = 0.0
    
    for step in 1:n_steps
        step_type = (step % 2 == 1) ? "W0" : "W1"
        
        if step <= 3 || step % 8 == 0
            println("  Passo $step: $step_type")
        end
        
        push!(step_markers, current_time)
        push!(step_labels, step_type)
        
        # Calcular parâmetros variáveis para este passo
        J_current = J_profile_func(current_time, J_base; J_kwargs...)
        τ_factor_current = τ_profile_func(current_time, τ_base_factor; τ_kwargs...)
        
        J_angular = 2π * J_current
        τ_current = π / (τ_factor_current * J_angular)
        dt = τ_current / (n_points_per_step - 1)
        
        # Criar Hamiltoniano para este passo
        H_step = create_step_hamiltonian(N, step_type, J_angular)
        
        # Evolução dentro do passo com parâmetros atualizados dinamicamente
        for i in 1:n_points_per_step
            # Recalcular parâmetros para cada sub-passo (permite variação mais suave)
            t_substep = current_time + (i-1) * dt
            J_substep = J_profile_func(t_substep, J_base; J_kwargs...)
            τ_factor_substep = τ_profile_func(t_substep, τ_base_factor; τ_kwargs...)
            
            # Armazenar evolução dos parâmetros
            push!(J_evolution, J_substep)
            push!(τ_evolution, τ_factor_substep)
            
            # Calcular probabilidades
            probs = calculate_all_probabilities(current_state, N)
            push!(time_points, t_substep)
            push!(prob_evolution, copy(probs))
            
            # Evoluir estado (exceto no último ponto)
            if i < n_points_per_step
                # Usar parâmetros do sub-passo atual para evolução
                J_angular_substep = 2π * J_substep
                H_substep = create_step_hamiltonian(N, step_type, J_angular_substep)
                U_step = time_evolve(H_substep, dt)
                current_state = apply!(current_state, U_step)
            end
        end
        
        current_time += τ_current
        
        # Debug periódico
        if step % 8 == 0
            final_probs = prob_evolution[end]
            max_prob = maximum(final_probs)
            max_qubit = argmax(final_probs)
            J_final = J_evolution[end]/1e6
            τ_final = τ_evolution[end]
            println("    → t=$(round(current_time*1e6, digits=1))μs: Máx P$max_qubit = $(round(max_prob, digits=3))")
            println("      J = $(round(J_final, digits=2)) MHz, τ_factor = $(round(τ_final, digits=1))")
        end
    end
    
    return (
        time_points = time_points,
        prob_evolution = prob_evolution,
        J_evolution = J_evolution,
        τ_evolution = τ_evolution,
        step_markers = step_markers,
        step_labels = step_labels
    )
end

# --- 5. FUNÇÃO DE ANÁLISE ---
function calculate_all_probabilities(state, N)
    """Calcula probabilidade de cada qubit estar em |1⟩"""
    probs = zeros(Float64, N)
    for i in 1:N
        z_op = put(N, i => Z)
        z_exp = real(expect(z_op, state))
        probs[i] = (1 - z_exp) / 2
    end
    return probs
end

function analyze_variable_results(simulation_results)
    """Analisa resultados da simulação com parâmetros variáveis"""
    time_points = simulation_results.time_points
    prob_evolution = simulation_results.prob_evolution
    J_evolution = simulation_results.J_evolution
    τ_evolution = simulation_results.τ_evolution
    
    prob_matrix = hcat(prob_evolution...)
    time_μs = time_points .* 1e6
    J_MHz = J_evolution ./ 1e6
    
    # Métricas
    qubits_activated = sum([maximum(prob_matrix[i, :]) > 0.1 for i in 1:N])
    initial_center = sum([i * prob_matrix[i, 1] for i in 1:N])
    final_center = sum([i * prob_matrix[i, end] for i in 1:N])
    propagation_distance = final_center - initial_center
    total_probs = sum(prob_matrix, dims=1)[:]
    avg_conservation = sum(total_probs) / length(total_probs)
    
    return (
        prob_matrix = prob_matrix,
        time_μs = time_μs,
        J_MHz = J_MHz,
        τ_evolution = τ_evolution,
        qubits_activated = qubits_activated,
        propagation_distance = propagation_distance,
        avg_conservation = avg_conservation,
        total_probs = total_probs,
        step_markers = simulation_results.step_markers,
        step_labels = simulation_results.step_labels
    )
end

# --- 6. CONFIGURAÇÕES DE EXEMPLO ---

println("\n" * "="^60)
println("🎯 CONFIGURAÇÕES DE EXEMPLO DISPONÍVEIS:")
println("="^60)

# Exemplo 1: J constante, τ com rampa linear
example_1 = (
    name = "J constante + τ rampa",
    J_func = constant_profile,
    τ_func = linear_ramp,
    J_kwargs = NamedTuple(),
    τ_kwargs = (t_max=1.0e-6, final_factor=3.0),
    description = "J fixo, τ cresce linearmente (pulsos ficam mais longos)"
)

# Exemplo 2: J com modulação senoidal, τ constante  
example_2 = (
    name = "J senoidal + τ constante",
    J_func = sinusoidal_modulation,
    τ_func = constant_profile,
    J_kwargs = (freq=3e6, amplitude=0.4),
    τ_kwargs = NamedTuple(),
    description = "J oscila senoidalmente, τ fixo (força variável)"
)

# Exemplo 3: J com degrau, τ com decaimento
example_3 = (
    name = "J degrau + τ decaimento",
    J_func = step_function,
    τ_func = exponential_decay,
    J_kwargs = (t_switch=6e-7, factor_before=1.0, factor_after=2.5),
    τ_kwargs = (decay_rate=2e6,),
    description = "J salta em t_switch, τ decai exponencialmente"
)

# Exemplo 4: J com pulso gaussiano, τ com chirp
example_4 = (
    name = "J gaussiano + τ chirp",
    J_func = gaussian_pulse,
    τ_func = chirp_modulation,
    J_kwargs = (t_center=4e-7, width=1.5e-7, amplitude=1.5),
    τ_kwargs = (freq_start=2e6, freq_end=8e6, t_max=1.2e-6),
    description = "J tem pico gaussiano, τ com frequência variável"
)

examples = [example_1, example_2, example_3, example_4]

for (i, ex) in enumerate(examples)
    println("$i. $(ex.name)")
    println("   $(ex.description)")
    println()
end

# --- 7. FUNÇÃO PARA EXECUTAR EXEMPLO ---
function run_example(example_num::Int; n_steps=24)
    """Executa um dos exemplos predefinidos"""
    
    if example_num < 1 || example_num > length(examples)
        error("Exemplo deve estar entre 1 e $(length(examples))")
    end
    
    ex = examples[example_num]
    println("\n🎯 EXECUTANDO EXEMPLO $example_num: $(ex.name)")
    println("📝 $(ex.description)")
    
    # Executar simulação
    results = simulate_variable_quantum_walk(
        ex.J_func, ex.τ_func,
        1.0e6, 10.0, n_steps, 25;
        J_kwargs=ex.J_kwargs, τ_kwargs=ex.τ_kwargs
    )
    
    # Analisar resultados
    analysis = analyze_variable_results(results)
    
    println("\n📊 RESULTADOS:")
    println("  Qubits ativados: $(analysis.qubits_activated) de $N")
    println("  Propagação: $(round(analysis.propagation_distance, digits=2)) qubits")
    println("  Conservação: $(round(analysis.avg_conservation, digits=4))")
    
    return results, analysis
end

# --- 8. VISUALIZAÇÃO ---
function plot_variable_results(results, analysis, title_suffix="")
    """Cria visualização completa dos resultados com parâmetros variáveis"""
    
    fig = Figure(size = (1400, 1000))
    
    # Heatmap principal
    ax_heatmap = Axis(fig[1, 1:3],
        title = "Caminhada Quântica - Parâmetros Variáveis $title_suffix",
        xlabel = "Tempo (μs)", ylabel = "Qubit",
        yticks = (1:N, ["Q$i" for i in 1:N]))

    prob_matrix_transposed = transpose(analysis.prob_matrix)
    
    hm = Makie.heatmap!(ax_heatmap, analysis.time_μs, 1:N, prob_matrix_transposed,
                  colormap = :plasma, colorrange = (0, 1))
    
    Colorbar(fig[1, 4], hm, label = "P(|1⟩)")
    
    # Marcar transições W0/W1
    step_times_μs = analysis.step_markers[2:end] .* 1e6
    for (i, t_marker) in enumerate(step_times_μs)
        if t_marker <= maximum(analysis.time_μs)
            vlines!(ax_heatmap, [t_marker], color=:white, linewidth=0.6, alpha=0.5)
        end
    end
    
    # Gráficos de evolução dos parâmetros
    ax_J = Axis(fig[2, 1], title = "Evolução J(t)", xlabel = "Tempo (μs)", ylabel = "J (MHz)")
    lines!(ax_J, analysis.time_μs, analysis.J_MHz, color = :red, linewidth=2)
    
    ax_tau = Axis(fig[2, 2], title = "Evolução τ_factor(t)", xlabel = "Tempo (μs)", ylabel = "τ factor")
    lines!(ax_tau, analysis.time_μs, analysis.τ_evolution, color = :blue, linewidth=2)
    
    ax_conservation = Axis(fig[2, 3], title = "Conservação", xlabel = "Tempo (μs)", ylabel = "Σ P")
    lines!(ax_conservation, analysis.time_μs, analysis.total_probs, color = :purple, linewidth=2)
    hlines!(ax_conservation, [1.0], color=:black, linestyle=:dash)
    
    # Evolução de qubits selecionados
    ax_qubits = Axis(fig[3, 1:3], title = "Evolução Probabilidades (Qubits Chave)", 
                     xlabel = "Tempo (μs)", ylabel = "P(|1⟩)")
    
    key_qubits = [1, 3, 5, 8]
    colors = [:red, :blue, :green, :purple]
    
    for (idx, qubit) in enumerate(key_qubits)
        lines!(ax_qubits, analysis.time_μs, analysis.prob_matrix[qubit, :],
               linewidth=2, color=colors[idx], label="Q$qubit")
    end
    axislegend(ax_qubits, position=:rt)
    
    return fig
end

# --- 9. EXECUÇÃO AUTOMÁTICA ---
println("🚀 Executando exemplo demonstrativo...")
println("   (Você pode executar outros com: run_example(1), run_example(2), etc.)")

# Executar exemplo 2 como demonstração
demo_results, demo_analysis = run_example(2, n_steps=24)

# Criar visualização
fig = plot_variable_results(demo_results, demo_analysis, "- Exemplo 2")
display(fig)

println("\n" * "="^60)
println("✨ COMO USAR:")
println("="^60)
println("• run_example(1): J constante + τ rampa")
println("• run_example(2): J senoidal + τ constante") 
println("• run_example(3): J degrau + τ decaimento")
println("• run_example(4): J gaussiano + τ chirp")
println("\n• Para criar seus próprios perfis, defina funções como:")
println("  function meu_perfil(t, base_value; param1=1.0)")
println("      return base_value * (1.0 + param1 * função_de_t)")
println("  end")
println("\n• Execute com: simulate_variable_quantum_walk(meu_perfil, outro_perfil, ...)")
println("="^60)