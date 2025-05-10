import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import statsmodels.api as sm
# import warnings
# warnings.filterwarnings('ignore')

def load_data(file_path):
    df = pd.read_csv(file_path, index_col=0)
    return df

def explore_data(df):
    print("Dataset shape:", df.shape)
    print("\nFeature data types:")
    print(df.dtypes)
    print("\nCheck for missing values:")
    print(df.isnull().sum())
    
    # Check class balance for target variable
    print("\nTarget variable distribution:")
    target_dist = df['is_effective'].value_counts(normalize=True) * 100
    print(target_dist)
    
    return df

def analyze_feature_correlations(df):
    numeric_cols = df.select_dtypes(include=['int64', 'float64']).columns.tolist()    
    corr_matrix = df[numeric_cols].corr(method='spearman')
    
    #plt.figure(figsize=(5,5))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', linewidths=0.5, fmt='.2f', annot_kws={"size": 8})
    plt.title('Feature Correlation Matrix')
    plt.tight_layout()
    plt.savefig('feature_correlation_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Show correlations with target variable
    target_corr = corr_matrix['is_effective'].sort_values(ascending=False)
    print("Correlations with target variable 'is_effective':")
    print(target_corr)
    
    return corr_matrix

def prepare_data(df):
    cols_to_keep = ["genetic_association", "GWAS_evidence"]
    X = df[cols_to_keep]
    y = df['is_effective']
    
    cat_cols = X.select_dtypes(include=['object', 'category']).columns.tolist()
    num_cols = X.select_dtypes(include=['int64', 'float64']).columns.tolist()
    
    scaler = StandardScaler()
    if num_cols:
        X[num_cols] = scaler.fit_transform(X[num_cols])
    
    if cat_cols:
        encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore', drop='first')
        encoded_cats = encoder.fit_transform(X[cat_cols])
        feature_names = encoder.get_feature_names_out(cat_cols)
        
        encoded_df = pd.DataFrame(
            encoded_cats, 
            columns=feature_names,
            index=X.index
        )
        
        X = X.drop(cat_cols, axis=1)
        X = pd.concat([X, encoded_df], axis=1)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.5, random_state=42, stratify=y
    )
    
    print(f"Training set: {X_train.shape}")
    print(f"Testing set: {X_test.shape}")
    
    return X_train, X_test, y_train, y_test

def evaluate_model(model, X_test, y_test):
    X_test_const = sm.add_constant(X_test)    
    y_prob = model.predict(X_test_const)
    y_pred = (y_prob > 0.5).astype(int)    
    
    # Calculate ROC AUC for model
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    roc_auc = auc(fpr, tpr)
    print(f"\nROC AUC: {roc_auc:.4f}")
    
    # Calculate baseline ROC AUC (random classifier)
    no_skill_prob = np.ones(len(y_test)) * np.mean(y_test)
    fpr_baseline, tpr_baseline, _ = roc_curve(y_test, no_skill_prob)
    roc_auc_baseline = auc(fpr_baseline, tpr_baseline)
    print(f"Baseline ROC AUC: {roc_auc_baseline:.4f}")
    
    # Calculate PR AUC for model
    precision, recall, _ = precision_recall_curve(y_test, y_prob)
    pr_auc = auc(recall, precision)
    print(f"PR AUC: {pr_auc:.4f}")
    
    # Calculate baseline PR AUC
    no_skill = len(y_test[y_test == 1]) / len(y_test)
    precision_baseline = np.ones_like(recall) * no_skill
    pr_auc_baseline = auc(recall, precision_baseline)
    print(f"Baseline PR AUC: {pr_auc_baseline:.4f}")
    
    # ROC curve with baseline comparison
    plt.figure(figsize=(5,5))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'Model ROC (AUC = {roc_auc:.2f})')
    plt.plot(fpr_baseline, tpr_baseline, color='green', lw=2, linestyle='--', 
             label=f'Baseline ROC (AUC = {roc_auc_baseline:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('roc_curve.png')
    
    # Plot Precision-Recall curve with baseline
    plt.figure(figsize=(5,5))
    plt.plot(recall, precision, color='blue', lw=2, 
             label=f'Model (AUC = {pr_auc:.2f})')
    plt.plot([0, 1], [no_skill, no_skill], linestyle='--', color='red', 
             label=f'No Skill ({no_skill:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend()
    plt.tight_layout()
    plt.savefig('precision_recall_curve.png')
 
    return y_pred, y_prob

def train_linear_model(X_train, y_train, alpha=0.01):
    X_train_const = sm.add_constant(X_train)    
    final_logit = sm.Logit(y_train, X_train_const)
    # Train model with L1 regularization (Lasso) and robust standard errors
    final_model = final_logit.fit_regularized(
        method='l1',
        alpha=alpha,
        cov_type='HC3',
        trim_mode='auto',
        maxiter=10_000
    )
    # Try L2 regularization (Ridge) if L1 fails to converge
    if not hasattr(final_model, 'params'):
        print("L1 regularization did not converge, trying L2...")
        final_model = final_logit.fit_regularized(
            method='l2',
            alpha=alpha,
            cov_type='HC3',
            trim_mode='auto',
            maxiter=10_000
        )
    train_pred_proba = final_model.predict(X_train_const)
    train_precision, train_recall, _ = precision_recall_curve(y_train, train_pred_proba)
    train_pr_auc = auc(train_recall, train_precision)
    print(f"Model training Precision-Recall AUC: {train_pr_auc:.4f}")
    fpr, tpr, _ = roc_curve(y_train, train_pred_proba)
    train_roc_auc = auc(fpr, tpr)
    print(f"Model training ROC AUC: {train_roc_auc:.4f}")
    print(f"Regularization: alpha={alpha}")
    return final_model

def analyze_coefficients(model, feature_names):
    feature_names = ['const'] + feature_names.tolist()
    coefficients = model.params    
    odds_ratios = np.exp(coefficients)
    coef_df = pd.DataFrame({
        'Feature': feature_names,
        'Coefficient': coefficients,
        'Odds_Ratio': odds_ratios,
        'HC_Std_Error': model.bse,
        'z_value': model.tvalues,
        'p_value': model.pvalues,
        'CI_Lower_95': model.conf_int(alpha=0.05).iloc[:, 0],
        'CI_Upper_95': model.conf_int(alpha=0.05).iloc[:, 1]
    })
    
    # Sort by absolute coefficient value (excluding constant)
    coef_df['Abs_Coefficient'] = np.abs(coef_df['Coefficient'])
    coef_df = coef_df.sort_values('Abs_Coefficient', ascending=False).drop('Abs_Coefficient', axis=1)
    
    # Add significance indicator
    coef_df['Significant'] = coef_df['p_value'] < 0.05
    
    # Print top significant coefficients (excluding constant)
    print("\nTop significant features (p < 0.05):")
    sig_features = coef_df[coef_df['Significant'] & (coef_df['Feature'] != 'const')]
    if len(sig_features) > 0:
        print(sig_features.head(10)[['Feature', 'Coefficient', 'p_value', 'Odds_Ratio', 'CI_Lower_95', 'CI_Upper_95']])
    else:
        print("No statistically significant features found.")
 
    # Count significant coefficients (excluding constant)
    n_significant = sum(coef_df['Significant'] & (coef_df['Feature'] != 'const'))
    n_features = sum(coef_df['Feature'] != 'const')
    print(f"\nNumber of significant coefficients (p < 0.05): {n_significant} out of {n_features}")
    
    # Plot top coefficients with confidence intervals (excluding constant)
    plt.figure(figsize=(5,5))
    # Get top 15 features by absolute coefficient value (excluding constant)
    top_coefs = coef_df[coef_df['Feature'] != 'const'].head(15)
    # Plot coefficients with confidence intervals
    colors = ['blue' if sig else 'gray' for sig in top_coefs['Significant']]
    for i, (coef, ci_lower, ci_upper, color) in enumerate(zip(
            top_coefs['Coefficient'], 
            top_coefs['CI_Lower_95'], 
            top_coefs['CI_Upper_95'],
            colors)):
        plt.errorbar(
            x=coef,
            y=i,
            xerr=[[coef - ci_lower], [ci_upper - coef]],
            fmt='o',
            capsize=5,
            color=color
        )
    # Add feature names as y-tick labels
    plt.yticks(range(len(top_coefs)), top_coefs['Feature'])
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    plt.title('Top Logistic Regression Coefficients with 95% CIs')
    plt.xlabel('Coefficient Value')
    plt.tight_layout()
    plt.savefig('top_coefficients.png', dpi=300, bbox_inches='tight')
    
    return coef_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_path', type=str, help='Path to the input CSV file')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save results')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    df = load_data(args.file_path)
    df = explore_data(df)    
    corr_matrix = analyze_feature_correlations(df)    
    X_train, X_test, y_train, y_test = prepare_data(df)    
    lr_model = train_linear_model(X_train, y_train)    
    lr_pred, lr_prob = evaluate_model(lr_model, X_test, y_test)    
    coef_df = analyze_coefficients(lr_model, X_train.columns)    
    
    results = {
        'logistic_regression': lr_model,
        'lr_coefficients': coef_df,
    }
    
    return results, df

if __name__ == "__main__":
    results, df = main()